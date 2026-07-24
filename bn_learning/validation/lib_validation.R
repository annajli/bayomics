#!/usr/bin/env Rscript
# =============================================================================
# lib_validation.R  —  shared machinery for the held-out flu-response validation
# -----------------------------------------------------------------------------
# The L5-only pipeline (scripts 01-04) is generalized here so EVERY learned
# network (clinical L1a/L1b, olink L2, clinical+olink L3a/L3b, whole-blood L4,
# cell-freq L5, pseudobulk L6, and the joint L_all) can be validated against the
# same held-out influenza-vaccine response.
#
# Three problems are solved once, network-agnostically:
#
#   (1) SUBJECT-ID RECOVERY.  Every fitted `<net>_data.rds` is an N x p frame
#       with NO subject identifier (rownames are 1..N). We recover identity
#       WITHOUT the missing CMV-backfill file: each continuous node maps (via the
#       network's node-mapping CSV) to a column of data/processed/bn_ready_baseline
#       .csv, which carries subject.subjectGuid. Those columns are subject-LOCAL
#       (raw labs / Olink NPX / gene-set scores / ALR ratios), i.e. NOT sample-
#       standardized, so a fitted row's values EXACTLY equal its subject's values.
#       We keep every node whose values provably align (raw, log, or log1p) and
#       exact-match each fitted row to its baseline subject on that vector.
#
#   (2) HELD-OUT TARGET ASSEMBLY.  flu_response_validation.csv is reduced to the
#       per-subject summaries used throughout: igg_seroconvert_any / _n,
#       igg_log2fc_mean, hai_peak_max (+ all per-antigen columns kept).
#
#   (3) STRUCTURAL (mediation) VALIDATION — the part that actually discriminates
#       between networks. Attaching the outcome as a leaf of the 3 roots is
#       network-INVARIANT (the internal DAG is irrelevant once the leaf's only
#       parents are the roots), so it is computed once as a cohort benchmark.
#       What differs per network is which root -> feature edges it LEARNED (from
#       baseline data alone, never seeing flu). Composing each learned root->cell
#       edge (sign A) with the independently-measured held-out feature->flu slope
#       (sign B) should reproduce the literature root->response direction. The
#       fraction that do — especially among features that are themselves flu-
#       relevant — is that network's structural-validation score.
#
# Everything is host-agnostic: paths resolve relative to bn_learning/ whether a
# caller runs from bn_learning/ or bn_learning/validation/.
# =============================================================================

suppressMessages({
  library(bnlearn)
  library(readr)
  library(dplyr)
})

# --- Paths -------------------------------------------------------------------
val_here_bn <- if (file.exists("config")) "." else ".."
VAL_MODELS  <- file.path(val_here_bn, "output", "models")
VAL_MAPS    <- file.path(val_here_bn, "output", "mappings")
VAL_OUT     <- file.path(val_here_bn, "output", "validation")
VAL_BASE    <- file.path(val_here_bn, "..", "data", "processed", "bn_ready_baseline.csv")
VAL_FLU     <- file.path(val_here_bn, "..", "data", "processed", "flu_response_validation.csv")
dir.create(VAL_OUT, showWarnings = FALSE, recursive = TRUE)

ROOTS <- c("age_group", "sex", "cmv")

# Canonical (data, models, node-mapping) triple per network. Uses the version
# each team pipeline last committed as its fitted artifact.
NETWORKS <- list(
  L1a  = list(id = "L1a_clinical_full_v3",      map = "L1a_clinical_full_v3_node_mapping.csv",
              label = "Clinical (full)",        modality = "clinical"),
  L1b  = list(id = "L1b_clinical_curated",      map = "L1b_clinical_curated_node_mapping.csv",
              label = "Clinical (curated)",     modality = "clinical"),
  L2   = list(id = "L2_olink_v2",               map = "L2_olink_v2_node_mapping.csv",
              label = "Olink (plasma protein)", modality = "olink"),
  L3a  = list(id = "L3a_clinical_full_olink",   map = "L3a_clinical_full_olink_node_mapping.csv",
              label = "Clinical(full)+Olink",   modality = "clinical+olink"),
  L3b  = list(id = "L3b_clinical_curated_olink",map = "L3b_clinical_curated_olink_node_mapping.csv",
              label = "Clinical(curated)+Olink",modality = "clinical+olink"),
  L4   = list(id = "L4_wb_pathways",            map = "L4_wb_pathways_node_mapping.csv",
              label = "Whole-blood pathways",   modality = "wb"),
  L5   = list(id = "L5_freq_immunophenotype",   map = "L5_freq_immunophenotype_node_mapping.csv",
              label = "Immune cell frequency",  modality = "freq"),
  L6   = list(id = "L6_pb_signaling",           map = "L6_pb_signaling_node_mapping.csv",
              label = "Pseudobulk signaling",   modality = "pb"),
  L_all= list(id = "L_all",                     map = "L_all_node_mapping.csv",
              label = "Joint (all modalities)", modality = "all")
)

# --- Held-out flu summaries (built once, shared) -----------------------------
# Returns subject.subjectGuid + per-antigen targets + the four per-subject
# summaries used everywhere. Seroconvert columns arrive as "True"/"False".
load_flu_summaries <- function() {
  val <- suppressMessages(read_csv(VAL_FLU, show_col_types = FALSE))
  igg_fc     <- grep("^val\\.igg_fc_",         names(val), value = TRUE)
  igg_log2fc <- grep("^val\\.igg_log2fc_",     names(val), value = TRUE)
  igg_sc     <- grep("^val\\.igg_seroconvert_", names(val), value = TRUE)
  hai_peak   <- grep("^val\\.hai_peak_",       names(val), value = TRUE)
  val[igg_sc] <- lapply(val[igg_sc], function(x) toupper(as.character(x)) == "TRUE")
  val <- val %>%
    mutate(
      igg_seroconvert_any = rowSums(across(all_of(igg_sc)), na.rm = TRUE) > 0,
      igg_seroconvert_n   = rowSums(across(all_of(igg_sc)), na.rm = TRUE),
      igg_log2fc_mean     = rowMeans(across(all_of(igg_log2fc)), na.rm = TRUE),
      hai_peak_max        = apply(across(all_of(hai_peak)), 1, max, na.rm = TRUE)
    )
  keep <- c("subject.subjectGuid", igg_fc, igg_log2fc, igg_sc, hai_peak,
            "igg_seroconvert_any", "igg_seroconvert_n",
            "igg_log2fc_mean", "hai_peak_max")
  val %>% select(all_of(keep))
}

# The four summary targets and whether higher = better response (used to keep
# the mediation direction bookkeeping consistent across targets).
FLU_TARGETS <- c("hai_peak_max", "igg_log2fc_mean")   # continuous, well-powered

# --- Node-mapping loader (tolerant of the two CSV schemas in the repo) --------
load_node_map <- function(map_file) {
  mp <- suppressMessages(read_csv(file.path(VAL_MAPS, map_file), show_col_types = FALSE))
  names(mp)[1:2] <- c("alias", "csv_column")
  mp$log_transformed <- if ("log_transformed" %in% names(mp))
    as.logical(mp$log_transformed) else FALSE
  mp$alr_transformed <- if ("alr_transformed" %in% names(mp))
    as.logical(mp$alr_transformed) else grepl("alr", mp$alias)
  mp
}

# For one candidate column, find the subject-local transform of the baseline
# column whose values reproduce the fitted values (raw / log / log1p), or NULL.
.detect_transform <- function(fitv, colv) {
  cand <- suppressWarnings(list(
    raw   = colv,
    log   = log(colv),
    log1p = log1p(colv)
  ))
  for (nm in names(cand)) {
    cv <- cand[[nm]]; cv <- cv[is.finite(cv)]
    if (length(cv) < length(fitv)) next
    ok <- vapply(fitv, function(x) any(abs(cv - x) < 1e-6 * (1 + abs(x))), logical(1))
    if (mean(ok) > 0.98) return(nm)
  }
  NULL
}

# =============================================================================
# recover_join(net) — recover subject IDs for a network and join flu outcomes.
# Returns the joined data.frame (also cached to <id>_validation_joined.rds/csv).
# =============================================================================
recover_join <- function(net, verbose = TRUE) {
  spec <- NETWORKS[[net]]
  d  <- readRDS(file.path(VAL_MODELS, paste0(spec$id, "_data.rds")))
  d  <- as.data.frame(d)
  mp <- load_node_map(spec$map)
  bn <- suppressMessages(read_csv(VAL_BASE, show_col_types = FALSE))

  cont <- setdiff(names(d), ROOTS)

  # Build the exact-match key from every continuous node that (a) maps to a
  # baseline column and (b) is not ALR (ALR ratios depend on a dropped, sample-
  # chosen reference so they need not appear verbatim) and (c) provably aligns
  # under a subject-local transform.
  key_fit <- list(); key_base <- list(); used <- character(0)
  for (nd in cont) {
    r <- which(mp$alias == nd); if (!length(r)) next
    csvc <- mp$csv_column[r[1]]
    if (isTRUE(mp$alr_transformed[r[1]])) next
    if (!(csvc %in% names(bn))) next
    colv <- as.numeric(bn[[csvc]]); fitv <- as.numeric(d[[nd]])
    tf <- .detect_transform(fitv, colv); if (is.null(tf)) next
    basev <- switch(tf, raw = colv, log = suppressWarnings(log(colv)), log1p = log1p(colv))
    key_fit[[nd]]  <- round(fitv, 6)
    key_base[[nd]] <- round(basev, 6)
    used <- c(used, nd)
  }
  if (length(used) < 3)
    stop(sprintf("[%s] only %d aligning columns — cannot recover subjects.",
                 net, length(used)))

  sig_fit  <- do.call(paste, c(key_fit,  sep = "|"))
  sig_base <- do.call(paste, c(key_base, sep = "|"))
  idx <- match(sig_fit, sig_base)
  if (anyNA(idx))
    stop(sprintf("[%s] %d/%d fitted rows unmatched.", net, sum(is.na(idx)), nrow(d)))
  if (anyDuplicated(idx))
    stop(sprintf("[%s] non-unique subject match — key not discriminating.", net))

  d$subject.subjectGuid <- bn$subject.subjectGuid[idx]

  # Cross-check the two roots we CAN verify against baseline (age, sex complete;
  # cmv was backfilled from a file that is missing, so it is not cross-checkable).
  base_roots <- bn %>%
    transmute(subject.subjectGuid,
              base_age = subject.ageGroup, base_sex = subject.biologicalSex)
  chk <- d %>% select(subject.subjectGuid, age_group, sex) %>%
    left_join(base_roots, by = "subject.subjectGuid")
  stopifnot(all(as.character(chk$age_group) == chk$base_age),
            all(as.character(chk$sex)       == chk$base_sex))

  flu <- load_flu_summaries()
  joined <- d %>% left_join(flu, by = "subject.subjectGuid")
  n_out <- sum(!is.na(joined$igg_seroconvert_any))

  if (verbose)
    cat(sprintf("[%-5s] recovered %d subjects on %d aligning cols; %d have flu outcomes\n",
                net, nrow(d), length(used), n_out))

  saveRDS(joined, file.path(VAL_OUT, paste0(spec$id, "_validation_joined.rds")))
  invisible(joined)
}

# =============================================================================
# leaf_root_query(joined) — network-INVARIANT cohort benchmark.
# Attach each held-out summary as a leaf of the 3 roots on an empty-internal DAG
# and read off the literature contrasts (Older-Young, Male-Female, CMV+ - CMV-).
# Identical for every network by construction, so run once on the largest join.
# =============================================================================
leaf_root_query <- function(joined) {
  d <- joined
  d$age_group <- factor(d$age_group, levels = c("Young Adult", "Older Adult"))
  d$sex       <- factor(d$sex,       levels = c("Female", "Male"))
  d$cmv       <- factor(d$cmv,       levels = c("Negative", "Positive"))
  contrasts_of <- function(ft) {
    f <- lm(reformulate(c("age_group", "sex", "cmv"), ft), data = d)
    co <- coef(f)
    data.frame(target = ft,
               older_minus_young = round(co[grep("age_group", names(co))[1]], 3),
               male_minus_female = round(co[grep("sex", names(co))[1]], 3),
               cmvpos_minus_neg  = round(co[grep("cmv", names(co))[1]], 3),
               row.names = NULL)
  }
  bind_rows(lapply(c("igg_seroconvert_any", "igg_log2fc_mean", "hai_peak_max"),
                   contrasts_of))
}

# =============================================================================
# structural_validation(net, joined) — the per-network discriminating test.
# For every confident root -> feature edge the network LEARNED (strength >= thr),
# compose sign(root->feature) x sign(feature->held-out-flu) and compare to the
# literature root->response direction (age Older -> lower, cmv Positive -> lower,
# sex Male -> lower). Returns a per-edge table plus a one-row network scorecard.
# =============================================================================
structural_validation <- function(net, joined, str_thr = 0.5, verbose = TRUE) {
  spec <- NETWORKS[[net]]
  m  <- readRDS(file.path(VAL_MODELS, paste0(spec$id, "_models.rds")))
  bs <- m$boot_str
  d  <- joined

  # Confident root -> feature edges (feature must be a continuous node we have).
  fr <- bs[bs$from %in% ROOTS & bs$strength >= str_thr & bs$direction >= 0.5, ]
  fr <- fr[fr$to %in% names(d), , drop = FALSE]
  fr <- fr[order(-fr$strength), ]

  lit_dir <- function(root) if (root == "sex") "lower(male)" else "lower"  # nonref -> lower
  nonref  <- function(root) switch(root, age_group = "Older Adult", sex = "Male", cmv = "Positive")

  rows <- list()
  for (i in seq_len(nrow(fr))) {
    root <- fr$from[i]; cell <- fr$to[i]
    grp <- d[[root]]; nr <- nonref(root)
    a <- mean(d[[cell]][grp == nr], na.rm = TRUE) - mean(d[[cell]][grp != nr], na.rm = TRUE)
    for (ft in FLU_TARGETS) {
      dd <- d; dd$cz <- scale(dd[[cell]])[, 1]
      f <- lm(reformulate("cz", ft), data = dd)
      b  <- coef(f)["cz"]; p <- summary(f)$coefficients["cz", 4]
      implied <- ifelse(sign(a) * sign(b) < 0, "lower", "higher")
      rows[[length(rows) + 1]] <- data.frame(
        network = net, root = root, feature = cell, strength = round(fr$strength[i], 3),
        flu_target = ft, A_root_to_feature = round(a, 3),
        B_feature_to_flu = round(b, 3), B_p = round(p, 3),
        implied_root_to_flu = implied,
        literature = sub("\\(.*", "", lit_dir(root)),
        agrees = implied == sub("\\(.*", "", lit_dir(root)),
        flu_relevant = p < 0.10,
        row.names = NULL)
    }
  }
  edges <- bind_rows(rows)

  # Network scorecard: overall agreement, agreement among flu-relevant edges, and
  # a strength x |B| weighted directional vote per root (aggregates weak signal).
  if (nrow(edges) == 0) {
    score <- data.frame(network = net, label = spec$label, n_subjects = nrow(d),
                        n_root_edges = 0, n_flu_relevant = 0,
                        agree_all = NA, agree_relevant = NA,
                        wvote_age = NA, wvote_sex = NA, wvote_cmv = NA)
    return(list(edges = edges, score = score))
  }
  wvote <- function(root) {
    e <- edges[edges$root == root, ]
    if (!nrow(e)) return(NA_real_)
    w <- e$strength * abs(e$B_feature_to_flu)
    v <- sum(w * ifelse(e$agrees, 1, -1)) / sum(w)  # +1 all-agree .. -1 all-disagree
    round(v, 3)
  }
  rel <- edges[edges$flu_relevant, ]
  score <- data.frame(
    network = net, label = spec$label, n_subjects = nrow(d),
    n_root_edges  = length(unique(paste(edges$root, edges$feature))),
    n_flu_relevant = length(unique(paste(rel$root, rel$feature))),
    agree_all      = round(mean(edges$agrees), 3),
    agree_relevant = if (nrow(rel)) round(mean(rel$agrees), 3) else NA_real_,
    wvote_age = wvote("age_group"), wvote_sex = wvote("sex"), wvote_cmv = wvote("cmv"))

  if (verbose)
    cat(sprintf("[%-5s] %d root-edges (%d flu-relevant); agree_all=%.0f%% agree_relevant=%s\n",
                net, score$n_root_edges, score$n_flu_relevant, 100 * score$agree_all,
                ifelse(is.na(score$agree_relevant), "NA",
                       paste0(round(100 * score$agree_relevant), "%"))))
  list(edges = edges, score = score)
}

# =============================================================================
# attenuation(joined, root, feature, target) — mediation attenuation test.
# Does adjusting the root -> outcome effect for the network-identified feature
# move it toward null (evidence the feature lies on the root -> response path)?
# An adjusted effect that flips sign or explodes is an unstable collinear fit,
# not genuine mediation, so it is flagged rather than reported as attenuation.
# =============================================================================
attenuation <- function(joined, root, feature, target = "hai_peak_max") {
  if (!feature %in% names(joined)) return(NULL)
  d <- joined; d$cz <- scale(d[[feature]])[, 1]
  base <- lm(reformulate(root, target), data = d)
  adj  <- lm(reformulate(c(root, "cz"), target), data = d)
  rt <- grep(root, names(coef(base)), value = TRUE)[1]
  eu <- coef(base)[rt]; ea <- coef(adj)[rt]
  data.frame(root = root, feature = feature, flu_target = target,
             effect_unadj = round(eu, 3), effect_adj = round(ea, 3),
             pct_attenuation = round(100 * (1 - ea / eu), 1),
             stable = sign(ea) == sign(eu) && abs(ea) <= 3 * abs(eu),
             row.names = NULL)
}
