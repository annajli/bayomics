# Held-out Validation Pipeline

Validates the learned Bayesian networks against the **held-out influenza
vaccine response** (IgG fold change / seroconversion / peak HAI), which was
deliberately excluded from structure learning. Backbone: the **L5**
cell-frequency network.

**Read `RESULTS.md` for the findings.** This file is the how-to.

## Prerequisites

- R with `bnlearn`, `dplyr`, `readr`, `tidyr`.
  - `bnlearn` needs R ≥ 4.4.1 for the current CRAN build. On R 4.4.0 install an
    older source build (e.g. `bnlearn_4.9.4`); if the compile fails with
    `library 'gfortran' not found`, point R at Homebrew's gfortran via
    `~/.R/Makevars`, e.g.
    `FLIBS=-L/opt/homebrew/Cellar/gcc/<ver>/lib/gcc/current -lgfortran -lquadmath`.
- Fitted L5 model at `../output/models/L5_freq_immunophenotype_models.rds`
  and its data frame `..._data.rds` (produced by the `bn_learning` pipeline).
- `../../data/processed/{bn_ready_baseline.csv, flu_response_validation.csv}`.

## Run

From `bn_learning/`:

```r
Rscript validation/run_all.R          # all four steps in order
# or individually:
Rscript validation/01_assemble_validation_data.R
Rscript validation/02_empirical_benchmark.R
Rscript validation/03_bn_leaf_queries.R
Rscript validation/04_mediation_check.R
```

Scripts resolve paths relative to `bn_learning/`, so they also run unchanged
from inside `validation/`.

## What each step does

| Script | Purpose | Key output |
|--------|---------|------------|
| `01_assemble_validation_data.R` | Recover subject IDs for the 88 modeled rows (exact feature match — no dependence on the missing CMV-backfill file), join held-out flu targets + antigen summaries | `L5_validation_joined.{rds,csv}` |
| `02_empirical_benchmark.R` | Establish the empirical ground truth: crude + baseline-adjusted root effects, per-antigen breakdown, literature scorecard | `emp_*.csv` |
| `03_bn_leaf_queries.R` | Attach the held-out outcome as a leaf of the roots on the fixed L5 DAG, fit its CPD, query `P(seroconvert \| age, cmv)` etc. by likelihood weighting | `bn_leaf_cpqueries.csv`, `bn_leaf_contrasts.csv` |
| `04_mediation_check.R` | Structural test: compose learned `root → cell` edges with held-out `cell → flu` relationships; mediation attenuation | `mediation_edges.csv`, `mediation_attenuation.csv` |

## Design in one paragraph

The flu response is not a node in any network. To produce conditional
probabilities we hold the learned L5 DAG fixed and attach the outcome as a leaf
(step 3) — the literal `P(seroconversion | age, CMV)`. Because outcome-on-roots
just reproduces the cohort's (baseline-confounded) association, the decisive
test is step 4: the network's `root → immune-cell` edges were learned without
ever seeing flu response, so composing them with the independently observed
`cell → flu` relationship in held-out data — and recovering the literature
`root → response` direction — is the real evidence the structure captures
biology. See `RESULTS.md`.
