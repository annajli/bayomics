# 00_install_packages.R
# One-time setup — run this script before any .Rmd notebooks
# R 4.5.2 required

# CRAN packages
install.packages(c(
  "bnlearn",
  "ggplot2",
  "dplyr",
  "readr",
  "tidyr",
  "knitr",
  "rmarkdown",
  "ggraph",
  "igraph",
  "here",
  "BiocManager"
))

# Bioconductor (for Rgraphviz DAG plotting)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Rgraphviz")

# Verify installation
cat("Checking packages...\n")
for (pkg in c("bnlearn", "Rgraphviz", "ggplot2", "dplyr", "readr",
              "tidyr", "knitr", "rmarkdown", "ggraph", "igraph",
              "here", "BiocManager")) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  OK: %s (%s)\n", pkg, packageVersion(pkg)))
  } else {
    cat(sprintf("  MISSING: %s\n", pkg))
  }
}
