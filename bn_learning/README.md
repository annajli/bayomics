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
| L1a_v2 | `config/L1a_clinical_full_v2.R` | 3 + 44 | L1a with BMI/height excluded (19.6% missing each; n rises 69→81) |
| L1b | `config/L1b_clinical_curated.R` | 3 + 22 | Literature-curated immune-aging panel (BMI excluded) |
| L1b_v2 | `config/L1b_clinical_curated_v2.R` | 3 + 22 | Same as L1b (BMI/height already excluded; identical results) |
| L2 | `config/L2_olink.R` | 3 + 35 | Olink plasma proteomics |
| L2_v2 | `config/L2_olink_v2.R` | 3 + 34 | Olink with MMP9 excluded (high missingness; n rises 72→88) |
| L3-L5 | TBD | TBD | Cross-modality and combined networks |

## Directory layout

```
bn_learning/
├── 00_install_packages.R        # One-time R package setup
├── 01_data_prep.Rmd             # Parameterized: params$config
├── 02_structure_learning.Rmd    # Parameterized: params$config
├── 03_model_comparison.Rmd      # Parameterized: params$version ("v1" or "v2")
├── config/                      # Network-level configs (sourced by 01/02)
├── output/
│   ├── models/                  # .rds model + data files
│   ├── figures/                 # Edge strength histograms, etc.
│   └── mappings/                # *_node_mapping.csv (short alias <-> bn_ready column)
└── rendered/                    # All HTML outputs + figure dirs
    ├── L1a/                     # 01/02 renders for L1a (v1 + v2)
    ├── L1b/                     # 01/02 renders for L1b (v1 + v2)
    ├── L2/                      # 01/02 renders for L2 (v1 + v2)
    └── comparison/              # 03 renders (v1, v2)
```

## Notebooks

| Notebook | Params | Purpose |
|----------|--------|---------|
| `01_data_prep.Rmd` | `config` | Load data, backfill CMV, transform, subset, complete-case |
| `02_structure_learning.Rmd` | `config` | Run hc/tabu/mmhc + bootstrap model averaging |
| `03_model_comparison.Rmd` | `version` | Compare algorithms, cross-level analysis, sanity queries |

Notebooks 01 and 02 are parameterized via `params$config`. Notebook 03 is
parameterized via `params$version` (`"v1"` or `"v2"`). Render examples:

```r
# Data prep + structure learning for a specific network level
rmarkdown::render("01_data_prep.Rmd", params = list(config = "config/L1b_clinical_curated.R"))

# Model comparison for v2 model files
rmarkdown::render("03_model_comparison.Rmd", params = list(version = "v2"))
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

Each network level has a `*_node_mapping.csv` in `output/mappings/` mapping short
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
- **CMV backfilled**: only 9/92 measured at baseline; 88/92 have CMV from other visits (stable trait); 4 subjects (BR1049, BR1056, BR1058, BR1059) have no CMV at any visit and are dropped at complete-case
- **Log-transformed**: CRP, ESR, ALT, AST, triglycerides, glucose (clinical networks only)
- **Olink NPX**: already log2-scale, no additional transform
- **Complete-case**: rows with any NA dropped (varies by network level)
- **NLR excluded**: deterministic function of neutrophil + lymphocyte counts
