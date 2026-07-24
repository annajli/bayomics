#!/usr/bin/env Rscript
# Run the full held-out validation pipeline in order.
# Usage (from bn_learning/):  Rscript validation/run_all.R
here_bn <- if (file.exists("config")) "." else ".."
scripts <- c("01_assemble_validation_data.R",   # L5: recover subjects + join flu
             "02_empirical_benchmark.R",         # L5: empirical ground truth
             "03_bn_leaf_queries.R",             # L5: P(response | roots) leaf query
             "04_mediation_check.R",             # L5: structural mediation
             "05_validate_all_networks.R")       # ALL networks: recover + structural
for (s in scripts) {
  cat("\n==================== ", s, " ====================\n")
  source(file.path(here_bn, "validation", s), local = new.env(), echo = FALSE)
}
cat("\nAll validation steps complete. See output/validation/ and validation/RESULTS.md\n")
