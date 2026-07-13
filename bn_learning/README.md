# Bayesian Network Structure Learning — BPGbio Capstone

## Overview

CLG (Conditional Linear Gaussian) Bayesian networks learned from the Sound Life
Cohort baseline multi-omics data (n=92). Three discrete root nodes
(age group, biological sex, CMV serostatus) drive continuous clinical and
molecular features.

## Network levels

| Level | Config | Nodes | Description |
|-------|--------|-------|-------------|
| L1a | `config/L1a_clinical_full.R` | 3 + 46 | All clinical variables (minus 6 exclude + NLR) |
| L1b | `config/L1b_clinical_curated.R` | 3 + 23 | Literature-curated immune-aging panel |
| L2 | `config/L2_olink.R` | 3 + 35 | Olink plasma proteomics |
| L3-L5 | TBD | TBD | Cross-modality and combined networks |

## Notebooks

| Notebook | Purpose |
|----------|---------|
| `01_data_prep.Rmd` | Load data, backfill CMV, transform, subset, complete-case |
| `02_structure_learning.Rmd` | Run hc/tabu/mmhc + bootstrap model averaging |
| `03_model_comparison.Rmd` | Compare algorithms, cross-level analysis, sanity queries |

All notebooks are parameterized via `params$config`. Render with:

```r
rmarkdown::render("01_data_prep.Rmd", params = list(config = "config/L1b_clinical_curated.R"))
```

## For the validator

### Loading a fitted model

```r
models <- readRDS("output/models/L1b_clinical_curated_models.rds")
fit <- models$fit_tabu          # fitted CLG parameters on tabu DAG
avg <- models$avg_opt           # averaged network (optimal threshold)

# Run a conditional probability query
library(bnlearn)
cpquery(fit,
        event = (crp > 1.0),
        evidence = (age_group == "Older Adult"),
        n = 50000)
```

### Node name mapping

Each network level has a `*_node_mapping.csv` in `output/` mapping short
aliases (used in the BN) back to `bn_ready_baseline.csv` column names.

### Validation data

`../data/processed/flu_response_validation.csv` — join on `subject.subjectGuid`.
Key targets: `val.igg_seroconvert`, `val.igg_fc`, `val.hai_peak` (per antigen).

### Algorithms

- `bn_hc`: hill-climbing (greedy)
- `bn_tabu`: tabu search (primary — escapes local optima)
- `bn_mmhc`: hybrid (constraint skeleton + score-based orientation)
- `boot_str`: bootstrap edge strengths (1000 replicates, tabu-based)
- `avg_opt` / `avg_85` / `avg_50`: averaged networks at optimal / 0.85 / 0.50 thresholds

### Key design decisions

- **CLG**: discrete roots cannot have continuous parents (enforced by blacklist)
- **CMV backfilled**: only 9/92 measured at baseline; all 92 have CMV from other visits (stable trait)
- **Log-transformed**: CRP, ESR, ALT, AST, triglycerides, glucose (clinical networks only)
- **Olink NPX**: already log2-scale, no additional transform
- **Complete-case**: rows with any NA dropped (varies by network level)
- **NLR excluded**: deterministic function of neutrophil + lymphocyte counts
