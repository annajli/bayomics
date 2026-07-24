# =============================================================================
# _paths.R — shared, location-agnostic path anchors for the validation scripts.
# Sourced by 01–05. Walks up from the working directory to the repo root (the
# dir containing both data/processed and bn_learning), so every script runs the
# same whether invoked from the repo root or from validation/.
# =============================================================================
find_repo_root <- function() {
  d <- normalizePath(getwd())
  for (i in 1:6) {
    if (dir.exists(file.path(d, "data", "processed")) &&
        dir.exists(file.path(d, "bn_learning"))) return(d)
    d <- dirname(d)
  }
  stop("Could not locate repo root (expected data/processed and bn_learning/).")
}
ROOT      <- find_repo_root()
BL_MODELS <- file.path(ROOT, "bn_learning", "output", "models")   # bootstrap models
BL_CONFIG <- file.path(ROOT, "bn_learning", "config")
DATA_PROC <- file.path(ROOT, "data", "processed")
VAL_DIR   <- file.path(ROOT, "validation")
OUT_DIR   <- file.path(VAL_DIR, "output")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
