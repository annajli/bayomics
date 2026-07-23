#!/usr/bin/env Rscript
# =============================================================================
# 01_assemble_validation_data.R
# -----------------------------------------------------------------------------
# HELD-OUT (VALIDATION) BENCHMARK — influenza vaccine response
#
# Purpose
#   Build the subject-aligned table that links each fitted-network subject to
#   its held-out flu-vaccine outcomes (IgG fold change, seroconversion, HAI
#   peak). The flu response was HELD OUT of Bayesian-network structure learning
#   (serology is only measured at vaccine visits, which would inject structured
#   missingness into the baseline network), so it lives in a separate file and
#   must be joined back in here for validation.
#
# The subject-ID recovery problem
#   The fitted L5 model data (output/models/L5_freq_immunophenotype_data.rds)
#   is an 88 x 91 frame with NO subject identifier — its rownames are just
#   1..88. To attach per-subject flu outcomes we must recover which 88 subjects
#   those rows are.
#
#   We do this WITHOUT needing the (now-missing) raw CMV-backfill metadata:
#     * L5's continuous nodes are cell-frequency ALR features, which are exact
#       deterministic transforms of the freq.* columns in bn_ready_baseline.csv
#       (0% missingness) — independent of CMV.
#     * We reconstruct those ALR features for all 92 baseline subjects (carrying
#       subject.subjectGuid), then EXACT-MATCH each fitted-model row to its
#       baseline subject on the feature vector.
#     * The discrete roots (age_group, sex, cmv — including the CMV backfill)
#       are then taken from the authoritative fitted-model frame.
#
# Output
#   output/validation/L5_validation_joined.rds  (and .csv)
#     one row per fitted-model subject (88), with: subject.subjectGuid, roots,
#     the L5 ALR features, and all held-out flu-response targets + summaries.
# =============================================================================

suppressMessages({
  library(readr)
  library(dplyr)
})

# --- Paths (script is host-agnostic: resolve relative to bn_learning/) --------
# Works whether invoked from bn_learning/ or bn_learning/validation/.
here_bn <- if (file.exists("config")) "." else ".."
cfg_path   <- file.path(here_bn, "config", "L5_freq_immunophenotype.R")
base_path  <- file.path(here_bn, "..", "data", "processed", "bn_ready_baseline.csv")
val_path   <- file.path(here_bn, "..", "data", "processed", "flu_response_validation.csv")
model_path <- file.path(here_bn, "output", "models", "L5_freq_immunophenotype_data.rds")
out_dir    <- file.path(here_bn, "output", "validation")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

stopifnot(file.exists(cfg_path), file.exists(base_path),
          file.exists(val_path), file.exists(model_path))

# --- 1. Config: NODE_MAP, PARENT_MAP, ROOTS, LOG_TRANSFORM --------------------
source(cfg_path)
cat(sprintf("Config: %s (%d continuous nodes)\n", NETWORK_ID, length(CONTINUOUS_MAP)))

# --- 2. Reconstruct the L5 modeling features for ALL 92 baseline subjects -----
# Mirrors bn_learning/01_data_prep_L5.Rmd sections 4, 5, 5b exactly, but keeps
# subject.subjectGuid alongside so we can recover identity by feature matching.
bn_raw <- read_csv(base_path, show_col_types = FALSE)
stopifnot(nrow(bn_raw) == 92)

# Select + rename continuous freq columns to their short aliases (roots handled
# separately — we take roots from the fitted frame, not reconstructed here).
cont_sel <- CONTINUOUS_MAP  # named: alias = csv_column
missing_cols <- setdiff(unname(cont_sel), names(bn_raw))
stopifnot(length(missing_cols) == 0)

feat <- bn_raw %>%
  select(subject.subjectGuid, all_of(unname(cont_sel)))
names(feat)[-1] <- names(cont_sel)                       # csv col -> alias
feat <- feat %>% mutate(across(-subject.subjectGuid, as.numeric))

# LOG_TRANSFORM is empty for L5 (freq already log-scale) but honor it generally.
for (v in LOG_TRANSFORM) feat[[v]] <- log1p(feat[[v]])

# --- 2b. Additive log-ratio (ALR) transform (config Section 5b) ---------------
# For each sibling group under a shared parent, and for the implicit top-level
# group (L1 lineages), replace all-but-the-highest-mean sibling with
# (sibling - reference) [log-scale] and drop the reference. Deterministic on
# freq columns, so it reproduces the fitted features exactly.
stopifnot(exists("PARENT_MAP"))
child_groups <- split(names(PARENT_MAP), unlist(PARENT_MAP))
top_level <- setdiff(unique(unlist(PARENT_MAP)), names(PARENT_MAP))
top_level <- intersect(top_level, names(feat))
if (length(top_level) >= 2) child_groups[["__top_level__"]] <- top_level

for (grp in names(child_groups)) {
  children <- intersect(child_groups[[grp]], names(feat))
  if (length(children) < 2) next
  means <- sapply(children, function(c) mean(feat[[c]], na.rm = TRUE))
  ref <- names(which.max(means))
  for (child in setdiff(children, ref)) {
    feat[[paste0(child, "_alr")]] <- feat[[child]] - feat[[ref]]
    feat[[child]] <- NULL
  }
  feat[[ref]] <- NULL
}
cat(sprintf("Reconstructed ALR features: %d subjects x %d features\n",
            nrow(feat), ncol(feat) - 1))

# --- 3. Load fitted-model frame and recover subject IDs by exact matching -----
model_df <- readRDS(model_path)
feat_cols <- setdiff(names(model_df), c("age_group", "sex", "cmv"))
stopifnot(all(feat_cols %in% names(feat)))          # same feature namespace

# Build a rounded signature per row on a handful of features -> robust exact key
# (full-vector match; rounding guards against float printback differences).
sig <- function(df, cols) apply(round(as.matrix(df[, cols]), 8), 1,
                                 paste, collapse = "|")
model_sig <- sig(model_df, feat_cols)
feat_sig  <- sig(feat,     feat_cols)

match_idx <- match(model_sig, feat_sig)
if (any(is.na(match_idx))) {
  stop(sprintf("Failed to match %d fitted rows to baseline subjects.",
               sum(is.na(match_idx))))
}
if (anyDuplicated(match_idx)) stop("Non-unique subject match — alignment ambiguous.")

model_df$subject.subjectGuid <- feat$subject.subjectGuid[match_idx]

# --- 3b. Verify alignment: matched feature values equal to 1e-8 ---------------
max_abs_diff <- max(abs(as.matrix(model_df[, feat_cols]) -
                        as.matrix(feat[match_idx, feat_cols])))
cat(sprintf("Alignment check: %d/%d rows matched, max |feature diff| = %.2e\n",
            sum(!is.na(match_idx)), nrow(model_df), max_abs_diff))
stopifnot(max_abs_diff < 1e-6)

# Sanity cross-check against baseline roots we CAN verify (age, sex are complete
# in bn_ready_baseline; cmv cannot be — it was backfilled from a missing file).
base_roots <- bn_raw %>%
  select(subject.subjectGuid,
         base_age = subject.ageGroup, base_sex = subject.biologicalSex)
chk <- model_df %>%
  select(subject.subjectGuid, age_group, sex) %>%
  left_join(base_roots, by = "subject.subjectGuid")
stopifnot(all(as.character(chk$age_group) == chk$base_age))
stopifnot(all(as.character(chk$sex) == chk$base_sex))
cat("Root cross-check OK: fitted age_group/sex match bn_ready_baseline.\n")

# --- 4. Load held-out flu-response targets and build summaries ----------------
val <- read_csv(val_path, show_col_types = FALSE)

igg_fc_cols       <- grep("^val\\.igg_fc_",       names(val), value = TRUE)  # 7 antigens
igg_log2fc_cols   <- grep("^val\\.igg_log2fc_",   names(val), value = TRUE)
igg_seroconv_cols <- grep("^val\\.igg_seroconvert_", names(val), value = TRUE)
hai_peak_cols     <- grep("^val\\.hai_peak_",     names(val), value = TRUE)  # 9 antigens
cat(sprintf("Targets: %d IgG-FC, %d seroconvert, %d HAI-peak antigens\n",
            length(igg_fc_cols), length(igg_seroconv_cols), length(hai_peak_cols)))

# Seroconvert columns are "True"/"False" strings -> logical
val[igg_seroconv_cols] <- lapply(val[igg_seroconv_cols],
                                 function(x) toupper(as.character(x)) == "TRUE")

# Per-subject SUMMARY targets across antigens:
#   * igg_seroconvert_any  : converts (>=4x) to >= 1 measured antigen (binary)
#   * igg_seroconvert_n    : number of antigens converted (0-7)
#   * igg_log2fc_mean      : mean log2 fold-change across antigens (continuous)
#   * hai_peak_max         : max peak HAI inhibition across antigens (continuous)
val_sum <- val %>%
  mutate(
    igg_seroconvert_any = rowSums(across(all_of(igg_seroconv_cols)), na.rm = TRUE) > 0,
    igg_seroconvert_n   = rowSums(across(all_of(igg_seroconv_cols)), na.rm = TRUE),
    igg_log2fc_mean     = rowMeans(across(all_of(igg_log2fc_cols)),  na.rm = TRUE),
    hai_peak_max        = apply(across(all_of(hai_peak_cols)), 1, max, na.rm = TRUE)
  )

# Keep IDs + all raw per-antigen targets + summaries
keep_cols <- c("subject.subjectGuid",
               igg_fc_cols, igg_log2fc_cols, igg_seroconv_cols, hai_peak_cols,
               "igg_seroconvert_any", "igg_seroconvert_n",
               "igg_log2fc_mean", "hai_peak_max")
val_keep <- val_sum %>% select(all_of(keep_cols))

# --- 5. Join outcomes onto the subject-aligned model frame --------------------
joined <- model_df %>%
  left_join(val_keep, by = "subject.subjectGuid")

n_no_outcome <- sum(is.na(joined$igg_seroconvert_any))
cat(sprintf("Joined outcomes: %d/%d subjects have flu-response data (%d missing)\n",
            sum(!is.na(joined$igg_seroconvert_any)), nrow(joined), n_no_outcome))

# --- 6. Save ------------------------------------------------------------------
saveRDS(joined, file.path(out_dir, "L5_validation_joined.rds"))
write_csv(joined, file.path(out_dir, "L5_validation_joined.csv"))
cat(sprintf("\nSaved: %s (%d x %d)\n",
            file.path(out_dir, "L5_validation_joined.rds"),
            nrow(joined), ncol(joined)))

# Quick orientation print
cat("\nRoot x seroconvert(any) cross-tab:\n")
print(table(age = joined$age_group, cmv = joined$cmv,
            seroconvert = joined$igg_seroconvert_any))
