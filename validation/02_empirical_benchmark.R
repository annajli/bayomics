#!/usr/bin/env Rscript
# =============================================================================
# 02_empirical_benchmark.R
# -----------------------------------------------------------------------------
# HELD-OUT VALIDATION — Step 2: empirical ground-truth benchmark
#
# Before asking whether the Bayesian network REPRODUCES literature-established
# flu-response relationships, we first establish what those relationships
# actually look like in THIS held-out cohort. The BN's inference is only
# meaningful as validation if we know the empirical direction/magnitude it is
# being checked against.
#
# Literature-established expectations (the directions the BN should recover):
#   * age_group = Older   -> SUPPRESSED response (immunosenescence): lower
#                            seroconversion, lower IgG fold-change, lower HAI.
#   * cmv = Positive       -> IMPAIRED response (CMV-driven memory inflation),
#                            classically strongest in older adults; modest/mixed.
#   * sex = Female         -> STRONGER antibody response than males.
#
# KNOWN CONFOUND — baseline titer. Fold-change / seroconversion are bounded by
# pre-vaccination (Day 0) titer: subjects who already sit high have little room
# to rise (a ceiling), which can INVERT the crude age effect (older adults,
# with more prior antigen exposure, differ at baseline). We therefore report
# BOTH crude stratified estimates AND baseline-adjusted regression estimates.
#
# Outputs (output/validation/):
#   emp_stratified_summary.csv   crude per-stratum means/rates
#   emp_regression_effects.csv   adjusted root effects (dir + p) per target
#   emp_per_antigen.csv          per-antigen age/cmv/sex effects
# =============================================================================

suppressMessages({ library(readr); library(dplyr); library(tidyr) })

source(if (file.exists("_paths.R")) "_paths.R" else "validation/_paths.R")
val_rds <- file.path(OUT_DIR, "L5_validation_joined.rds")
out_dir <- OUT_DIR
stopifnot(file.exists(val_rds))
d <- readRDS(val_rds)

# Ensure factor references match the network (Young/Female/Negative = baseline).
d$age_group <- factor(d$age_group, levels = c("Young Adult", "Older Adult"))
d$sex       <- factor(d$sex,       levels = c("Female", "Male"))
d$cmv       <- factor(d$cmv,       levels = c("Negative", "Positive"))

# Per-antigen column groups
igg_log2fc_cols   <- grep("^val\\.igg_log2fc_",     names(d), value = TRUE)
igg_seroconv_cols <- grep("^val\\.igg_seroconvert_", names(d), value = TRUE)
hai_peak_cols     <- grep("^val\\.hai_peak_",       names(d), value = TRUE)
igg_d0_cols       <- NULL  # baseline titers live in the raw file; see below

# Baseline (Day-0) IgG titers are needed as the confounder control. They are in
# the source validation csv, not carried into the joined frame — pull them.
val_raw <- read_csv(file.path(DATA_PROC, "flu_response_validation.csv"),
                    show_col_types = FALSE)
# The 7 real Day-0 titer columns share the antigen suffixes of the log2FC
# columns; the same "igg_d0" prefix also matches boolean QC flags
# (_replspread_flag, _hightiter_flag) which we must exclude.
antigens <- sub("^val\\.igg_log2fc_", "", igg_log2fc_cols)
d0_cols  <- paste0("val.igg_d0_", antigens)
stopifnot(all(d0_cols %in% names(val_raw)))
# mean log2 baseline titer across the 7 IgG antigens, per subject
val_raw$igg_d0_log2mean <- rowMeans(log2(val_raw[d0_cols]), na.rm = TRUE)
stopifnot(all(is.finite(val_raw$igg_d0_log2mean)))
d <- d %>% left_join(val_raw %>% select(subject.subjectGuid, igg_d0_log2mean),
                     by = "subject.subjectGuid")

cat(sprintf("Loaded %d subjects. Roots: age(%s) sex(%s) cmv(%s)\n\n",
            nrow(d),
            paste(table(d$age_group), collapse="/"),
            paste(table(d$sex), collapse="/"),
            paste(table(d$cmv), collapse="/")))

# =============================================================================
# 1. Crude stratified summary (one row per root level)
# =============================================================================
strat_one <- function(df, var) {
  df %>% group_by(level = .data[[var]]) %>%
    summarise(
      root = var,
      n = n(),
      seroconvert_any_rate = mean(igg_seroconvert_any),
      seroconvert_n_mean   = mean(igg_seroconvert_n),
      igg_log2fc_mean      = mean(igg_log2fc_mean, na.rm = TRUE),
      hai_peak_max_mean    = mean(hai_peak_max, na.rm = TRUE),
      igg_d0_log2mean      = mean(igg_d0_log2mean, na.rm = TRUE),
      .groups = "drop") %>%
    relocate(root)
}
strat <- bind_rows(strat_one(d, "age_group"),
                   strat_one(d, "sex"),
                   strat_one(d, "cmv"))
cat("=== Crude stratified summary ===\n")
print(as.data.frame(strat), digits = 3)
write_csv(strat, file.path(out_dir, "emp_stratified_summary.csv"))

# =============================================================================
# 2. Baseline-adjusted regression effects
#    seroconvert_any ~ age + sex + cmv           (logistic)
#    igg_log2fc_mean ~ age + sex + cmv (+ d0)     (linear, with/without baseline)
#    hai_peak_max    ~ age + sex + cmv            (linear)
#    Report each root's coefficient sign (direction) and p-value.
# =============================================================================
tidy_glm <- function(fit, target, model_label) {
  s <- summary(fit)$coefficients
  keep <- rownames(s) %in% c("age_groupOlder Adult", "sexMale", "cmvPositive")
  data.frame(
    target = target, model = model_label,
    term = rownames(s)[keep],
    estimate = s[keep, 1],
    p_value  = s[keep, ncol(s)],
    row.names = NULL)
}

eff <- bind_rows(
  tidy_glm(glm(igg_seroconvert_any ~ age_group + sex + cmv,
               family = binomial, data = d),
           "seroconvert_any", "logistic (crude)"),
  tidy_glm(glm(igg_seroconvert_any ~ age_group + sex + cmv + igg_d0_log2mean,
               family = binomial, data = d),
           "seroconvert_any", "logistic + baseline"),
  tidy_glm(lm(igg_log2fc_mean ~ age_group + sex + cmv, data = d),
           "igg_log2fc_mean", "linear (crude)"),
  tidy_glm(lm(igg_log2fc_mean ~ age_group + sex + cmv + igg_d0_log2mean, data = d),
           "igg_log2fc_mean", "linear + baseline"),
  tidy_glm(lm(hai_peak_max ~ age_group + sex + cmv, data = d),
           "hai_peak_max", "linear (crude)")
)
# Human-readable direction relative to the reference (Young / Female / Negative)
eff$direction <- ifelse(eff$estimate < 0, "lower", "higher")
cat("\n=== Baseline-adjusted root effects (relative to Young/Female/CMV-neg) ===\n")
print(eff, digits = 3)
write_csv(eff, file.path(out_dir, "emp_regression_effects.csv"))

# =============================================================================
# 3. Per-antigen effects (does the direction hold antigen-by-antigen?)
#    For each IgG antigen: age/cmv/sex effect on log2 FC + seroconvert rate.
# =============================================================================
antigen_of <- function(col, prefix) sub(prefix, "", col)
per_ant <- lapply(igg_log2fc_cols, function(col) {
  ant <- antigen_of(col, "val.igg_log2fc_")
  sc_col <- paste0("val.igg_seroconvert_", ant)
  dd <- d; dd$y <- dd[[col]]
  fit <- lm(y ~ age_group + sex + cmv, data = dd)
  co <- coef(fit)
  data.frame(
    antigen = ant,
    n_seroconvert = sum(dd[[sc_col]], na.rm = TRUE),
    mean_log2fc = mean(dd[[col]], na.rm = TRUE),
    age_older_effect = unname(co["age_groupOlder Adult"]),
    cmv_pos_effect   = unname(co["cmvPositive"]),
    sex_male_effect  = unname(co["sexMale"]))
})
per_ant <- bind_rows(per_ant)
cat("\n=== Per-antigen IgG log2FC: root effects (neg = lower than reference) ===\n")
print(per_ant, digits = 3)
write_csv(per_ant, file.path(out_dir, "emp_per_antigen.csv"))

# =============================================================================
# 4. Direction-vs-literature scorecard (printed)
# =============================================================================
cat("\n=== Direction vs. literature (baseline-adjusted models) ===\n")
lit <- data.frame(
  term = c("age_groupOlder Adult", "cmvPositive", "sexMale"),
  literature = c("lower (immunosenescence)",
                 "lower (CMV impairment)",
                 "lower (males weaker than females)"))
score <- eff %>%
  filter(model %in% c("logistic + baseline", "linear + baseline")) %>%
  left_join(lit, by = "term") %>%
  mutate(expected_dir = sub(" .*", "", literature),
         agrees = direction == expected_dir) %>%
  select(target, term, estimate, direction, expected_dir, agrees, p_value)
print(score, digits = 3)
write_csv(score, file.path(out_dir, "emp_literature_scorecard.csv"))

cat("\nDone. Wrote emp_*.csv to output/validation/\n")
