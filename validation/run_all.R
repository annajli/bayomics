#!/usr/bin/env Rscript
# Run the full held-out validation pipeline in order.
# Usage (from the repo root):  Rscript validation/run_all.R
#    (or from validation/):     Rscript run_all.R
val_dir <- if (file.exists("run_all.R")) "." else "validation"
scripts <- c("01_assemble_validation_data.R",   # L5: recover subjects + join flu
             "02_empirical_benchmark.R",         # L5: empirical ground truth
             "03_bn_leaf_queries.R",             # L5: P(response | roots) leaf query
             "04_mediation_check.R",             # L5: structural mediation
             "05_validate_all_networks.R")       # ALL networks x {bootstrap, cv}
for (s in scripts) {
  cat("\n==================== ", s, " ====================\n")
  source(file.path(val_dir, s), local = new.env(), echo = FALSE)
}
cat("\nAll validation steps complete. See validation/output/ and validation/RESULTS*.md\n")
