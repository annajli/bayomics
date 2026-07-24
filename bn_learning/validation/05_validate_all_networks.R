#!/usr/bin/env Rscript
# =============================================================================
# 05_validate_all_networks.R
# -----------------------------------------------------------------------------
# Held-out flu-response validation for EVERY learned network (not just L5).
#
# For each network in lib_validation::NETWORKS this:
#   1. recovers subject IDs and joins the held-out flu outcomes,
#   2. runs the structural (mediation) validation — the test that actually
#      distinguishes networks, since each learned different root -> feature edges,
#   3. contributes to a combined scorecard.
# The root-level leaf query is network-invariant, so it is computed once on the
# largest join as a cohort benchmark.
#
# Outputs (output/validation/):
#   <id>_validation_joined.rds          per-network subject-aligned join (cached)
#   all_networks_scorecard.csv          one row per network: agreement + w-votes
#   all_networks_structural_edges.csv   every root->feature edge, scored
#   all_networks_flagship_axes.csv      canonical marker axes, where present
#   leaf_root_benchmark.csv             network-invariant P/E(response | roots)
# =============================================================================

here_bn <- if (file.exists("config")) "." else ".."
source(file.path(here_bn, "validation", "lib_validation.R"))

nets <- names(NETWORKS)
cat("Validating", length(nets), "networks:", paste(nets, collapse = ", "), "\n\n")

all_edges <- list(); all_scores <- list(); joins <- list()
for (net in nets) {
  j <- tryCatch(recover_join(net), error = function(e) {
    cat(sprintf("[%-5s] SKIPPED: %s\n", net, conditionMessage(e))); NULL })
  if (is.null(j)) next
  joins[[net]] <- j
  sv <- structural_validation(net, j)
  all_edges[[net]]  <- sv$edges
  all_scores[[net]] <- sv$score
}

scorecard <- bind_rows(all_scores)
edges_all <- bind_rows(all_edges)
write_csv(scorecard, file.path(VAL_OUT, "all_networks_scorecard.csv"))
write_csv(edges_all, file.path(VAL_OUT, "all_networks_structural_edges.csv"))

# --- Flagship canonical axes (literature-anchored markers), where a network
#     learned the root->feature edge AND carries the feature. These are the most
#     interpretable single-edge validations. ------------------------------------
flagship <- tribble(
  ~root,        ~feature_pattern,                       ~marker,
  "age_group",  "^GDF15$",                              "GDF15 (canonical aging biomarker)",
  "age_group",  "^CCL11$",                              "CCL11/eotaxin (inflammaging)",
  "age_group",  "sox4_naive_cd4",                       "naive CD4 T (thymic output)",
  "sex",        "^LEP$",                                "leptin (sex-dimorphic)",
  "sex",        "creatinine",                           "creatinine (muscle mass)",
  "cmv",        "klrf1_gzmb_cd27_memory_cd4",           "cytotoxic memory CD4 (CMV inflation)",
  "cmv",        "adaptive_nk",                          "adaptive NK (CMV imprint)",
  "cmv",        "^KLRD1$",                              "KLRD1/NKG2D (CMV imprint)"
)
flag_rows <- list()
for (net in names(joins)) {
  e <- all_edges[[net]]; if (is.null(e) || !nrow(e)) next
  for (i in seq_len(nrow(flagship))) {
    hit <- e[e$root == flagship$root[i] &
             grepl(flagship$feature_pattern[i], e$feature) &
             e$flu_target == "hai_peak_max", ]
    if (nrow(hit)) flag_rows[[length(flag_rows) + 1]] <-
      cbind(marker = flagship$marker[i], hit)
  }
}
flagship_tbl <- if (length(flag_rows)) bind_rows(flag_rows) else data.frame()
write_csv(flagship_tbl, file.path(VAL_OUT, "all_networks_flagship_axes.csv"))

# --- Mediation attenuation for the three strongest continuous axes, in every
#     network that carries the feature. If adjusting the root->HAI effect for the
#     network-identified feature moves it toward null, the feature lies on the
#     root -> response path. --------------------------------------------------
atten_axes <- tribble(
  ~root,       ~feat_pat,
  "age_group", "sox4_naive_cd4",
  "age_group", "^GDF15$",
  "cmv",       "adaptive_nk"
)
atten_rows <- list()
for (net in names(joins)) {
  j <- joins[[net]]
  for (i in seq_len(nrow(atten_axes))) {
    feat <- grep(atten_axes$feat_pat[i], names(j), value = TRUE)
    feat <- setdiff(feat, ROOTS)
    if (!length(feat)) next
    a <- attenuation(j, atten_axes$root[i], feat[1])
    if (!is.null(a)) atten_rows[[length(atten_rows) + 1]] <- cbind(network = net, a)
  }
}
atten_tbl <- if (length(atten_rows)) bind_rows(atten_rows) else data.frame()
write_csv(atten_tbl, file.path(VAL_OUT, "all_networks_attenuation.csv"))

# --- Network-invariant leaf root benchmark (compute once on the biggest join) --
biggest <- joins[[ which.max(vapply(joins, nrow, integer(1))) ]]
leaf_bench <- leaf_root_query(biggest)
write_csv(leaf_bench, file.path(VAL_OUT, "leaf_root_benchmark.csv"))

# --- Console summary ---------------------------------------------------------
cat("\n===================== SCORECARD =====================\n")
print(scorecard, row.names = FALSE)
cat("\n============ FLAGSHIP CANONICAL AXES (HAI peak) ============\n")
if (nrow(flagship_tbl)) {
  print(flagship_tbl[, c("network","marker","A_root_to_feature",
                         "B_feature_to_flu","B_p","implied_root_to_flu",
                         "literature","agrees")], row.names = FALSE)
} else cat("(none matched)\n")
cat("\n============ MEDIATION ATTENUATION (root->HAI peak, adj. for feature) ============\n")
if (nrow(atten_tbl)) print(atten_tbl, row.names = FALSE) else cat("(none)\n")
cat("\n============ LEAF ROOT BENCHMARK (network-invariant) ============\n")
print(leaf_bench, row.names = FALSE)

cat("\nWrote: all_networks_scorecard.csv, all_networks_structural_edges.csv,\n",
    "       all_networks_flagship_axes.csv, leaf_root_benchmark.csv,\n",
    "       and <id>_validation_joined.rds for each network.\n")
