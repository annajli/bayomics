#!/usr/bin/env Rscript
# Run the full held-out validation pipeline in order.
# Usage (from bn_learning/):  Rscript validation/run_all.R
here_bn <- if (file.exists("config")) "." else ".."
scripts <- c("01_assemble_validation_data.R",
             "02_empirical_benchmark.R",
             "03_bn_leaf_queries.R",
             "04_mediation_check.R")
for (s in scripts) {
  cat("\n==================== ", s, " ====================\n")
  source(file.path(here_bn, "validation", s), local = new.env(), echo = FALSE)
}
cat("\nAll validation steps complete. See output/validation/ and validation/RESULTS.md\n")
