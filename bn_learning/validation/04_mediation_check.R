#!/usr/bin/env Rscript
# =============================================================================
# 04_mediation_check.R
# -----------------------------------------------------------------------------
# HELD-OUT VALIDATION — Step 4: structural (mediation) validation
#
# This is the strongest test of whether the network's LEARNED STRUCTURE captures
# real vaccine-response biology. The flu outcome was never in the network, yet
# the network learned root -> immune-cell edges purely from baseline data. If
# those edges encode real biology, then composing them with the (independently
# measured) cell -> flu-response relationship in the held-out data should
# reproduce the literature-established root -> response direction.
#
# For each strong root -> cell edge the L5 network learned:
#   (A) root -> cell : sign of the root's effect on the cell feature
#                      (network edge direction, confirmed against the data).
#   (B) cell -> flu  : sign of the cell feature's association with the held-out
#                      flu response (regression on the 88 joined subjects).
#   (C) implied root -> flu = sign(A) * sign(B), compared to literature.
#
# The biology we expect (mechanistic priors, independent of this dataset):
#   * age -> naive CD4 T cells (SOX4+): DOWN with age (thymic involution).
#     Naive CD4 help is required to mount new antibody responses -> naive UP
#     should mean BETTER response. Composed: age -> lower response.
#   * cmv -> cytotoxic memory CD4 (KLRF1/GZMB/CD27) and adaptive NK: UP with
#     CMV (memory inflation). Memory inflation crowds the repertoire and is
#     associated with WORSE vaccine response. Composed: cmv -> lower response.
#   * sex -> early memory B cells: the B-cell maturation axis underlying the
#     stronger female antibody response.
#
# Outputs (output/validation/):
#   mediation_edges.csv        per root->cell edge: A, B, implied, literature
#   mediation_attenuation.csv  root->outcome effect before/after adjusting for
#                              the network-identified mediating cell
# =============================================================================

suppressMessages({ library(bnlearn); library(dplyr); library(readr) })

here_bn <- if (file.exists("config")) "." else ".."
out_dir <- file.path(here_bn, "output", "validation")
m <- readRDS(file.path(here_bn, "output", "models",
                       "L5_freq_immunophenotype_models.rds"))
d <- readRDS(file.path(out_dir, "L5_validation_joined.rds"))

d$age_group <- factor(d$age_group, levels = c("Young Adult", "Older Adult"))
d$sex       <- factor(d$sex,       levels = c("Female", "Male"))
d$cmv       <- factor(d$cmv,       levels = c("Negative", "Positive"))

# The strong root->cell edges the network learned (from boot_str, strength>=0.5),
# each annotated with the mechanistic prior for cell -> response.
edges <- tribble(
  ~root,        ~cell,                                       ~boot_strength, ~cell_meaning,                              ~prior_cell_to_flu,
  "age_group",  "l3_sox4_naive_cd4_t_cell_alr",              0.816,          "naive CD4 T cells (thymic output)",        "higher = better (naive help)",
  "cmv",        "l3_klrf1_gzmb_cd27_memory_cd4_t_cell_alr",  1.000,          "cytotoxic memory CD4 (CMV inflation)",     "higher = worse (memory inflation)",
  "cmv",        "l3_adaptive_nk_cell_alr",                   0.800,          "adaptive/memory NK (CMV imprint)",         "higher = worse (CMV imprint)",
  "sex",        "l3_early_memory_b_cell_alr",                0.702,          "early memory B cells",                     "higher = better (B maturation)",
  "age_group",  "l2_cd16_monocyte_alr",                      0.661,          "CD16 non-classical monocytes",             "higher = worse (inflammaging)"
)

# Continuous flu targets used for the cell->flu regressions (binary is
# underpowered at n=88). hai_peak_max is the baseline-robust absolute measure;
# igg_log2fc_mean is the fold-change measure.
flu_targets <- c("hai_peak_max", "igg_log2fc_mean")

nonref_of <- function(root) switch(root,
  age_group = "Older Adult", sex = "Male", cmv = "Positive")

# ---------------------------------------------------------------------------
# (A) root -> cell : mean difference (non-reference minus reference level).
#     The network edge asserts a dependence; we report its data-confirmed sign.
# (B) cell -> flu  : standardized slope from lm(flu ~ scale(cell)).
# (C) implied and agreement vs literature (older/CMV+ -> lower; here we track
#     direction on the response scale).
# ---------------------------------------------------------------------------
rows <- list()
for (i in seq_len(nrow(edges))) {
  e <- edges[i, ]
  cell <- e$cell
  stopifnot(cell %in% names(d))
  nonref <- nonref_of(e$root)

  # (A) root effect on the cell
  grp <- d[[e$root]]
  a_diff <- mean(d[[cell]][grp == nonref], na.rm = TRUE) -
            mean(d[[cell]][grp != nonref], na.rm = TRUE)

  for (ft in flu_targets) {
    dd <- d
    dd$cell_z <- scale(dd[[cell]])[, 1]
    fit <- lm(reformulate("cell_z", ft), data = dd)
    b_slope <- coef(fit)["cell_z"]
    b_p     <- summary(fit)$coefficients["cell_z", 4]

    implied <- sign(a_diff) * sign(b_slope)   # implied root(nonref) -> flu sign
    rows[[length(rows) + 1]] <- data.frame(
      root = e$root, nonref_level = nonref, cell = cell,
      cell_meaning = e$cell_meaning, boot_strength = e$boot_strength,
      flu_target = ft,
      A_root_to_cell = round(a_diff, 3),
      B_cell_to_flu_slope = round(b_slope, 3),
      B_p = round(b_p, 3),
      prior_cell_to_flu = e$prior_cell_to_flu,
      implied_root_to_flu = ifelse(implied < 0, "lower", "higher"),
      literature_root_to_flu = ifelse(e$root == "sex", "lower(male)", "lower"),
      stringsAsFactors = FALSE)
  }
}
med <- bind_rows(rows)
med$agrees_literature <- med$implied_root_to_flu ==
  sub("\\(.*", "", med$literature_root_to_flu)

cat("=== Structural mediation: root -> cell (network) x cell -> flu (held-out) ===\n")
print(med[, c("root","cell_meaning","flu_target","A_root_to_cell",
              "B_cell_to_flu_slope","B_p","implied_root_to_flu",
              "literature_root_to_flu","agrees_literature")],
      row.names = FALSE)
write_csv(med, file.path(out_dir, "mediation_edges.csv"))

# ---------------------------------------------------------------------------
# Attenuation check: does adjusting for the network-identified mediating cell
# move the root -> outcome effect toward null (evidence the cell lies on the
# root -> response path)? Reported for the primary immunosenescence edge and
# the primary CMV edge, on hai_peak_max (the deconfounded target).
# ---------------------------------------------------------------------------
atten_one <- function(root, cell, ft) {
  dd <- d; dd$cell_z <- scale(dd[[cell]])[, 1]
  base <- lm(reformulate(root, ft), data = dd)
  adj  <- lm(reformulate(c(root, "cell_z"), ft), data = dd)
  rootterm <- grep(root, names(coef(base)), value = TRUE)[1]
  eu <- coef(base)[rootterm]; ea <- coef(adj)[rootterm]
  # An adjusted effect that flips sign or explodes in magnitude relative to the
  # unadjusted one is not a meaningful "attenuation" — it signals an unstable,
  # collinear/noisy mediator fit rather than genuine mediation. Flag it.
  stable <- sign(ea) == sign(eu) && abs(ea) <= 3 * abs(eu)
  data.frame(root = root, cell = cell, flu_target = ft,
             root_effect_unadj = round(eu, 3),
             root_effect_adj   = round(ea, 3),
             pct_attenuation = round(100 * (1 - ea / eu), 1),
             stable_estimate = stable,
             row.names = NULL)
}
atten <- bind_rows(
  atten_one("age_group", "l3_sox4_naive_cd4_t_cell_alr", "hai_peak_max"),
  atten_one("cmv", "l3_klrf1_gzmb_cd27_memory_cd4_t_cell_alr", "hai_peak_max"),
  atten_one("cmv", "l3_adaptive_nk_cell_alr", "hai_peak_max")
)
cat("\n=== Mediation attenuation of root->HAI-peak effect (adjusting for cell) ===\n")
print(atten, row.names = FALSE)
write_csv(atten, file.path(out_dir, "mediation_attenuation.csv"))

cat("\nDone. Wrote mediation_edges.csv and mediation_attenuation.csv\n")
