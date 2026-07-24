# Held-out Validation Pipeline

Top-level, cross-cutting validation of the learned Bayesian networks against the
**held-out influenza vaccine response** (IgG fold change / seroconversion / peak
HAI), which was deliberately excluded from structure learning. It covers **both**
structure-learning pipelines:

- **`bn_learning/`** — bootstrap model averaging (`source = "bootstrap"`)
- **`bn_cv/`** — cross-validated (`source = "cv"`)

Two layers:

- **L5 backbone** (scripts 01–04): the deep-dive on the cell-frequency network
  (bootstrap models). **Findings → `RESULTS.md`.**
- **All networks × both sources** (script 05 + `lib_validation.R`): the same
  validation applied to every network — clinical (L1a/L1b), Olink (L2),
  clinical+Olink (L3a/L3b), whole-blood (L4), cell freq (L5), pseudobulk (L6),
  and the joint **L_all** — under both bootstrap and CV.
  **Findings → `RESULTS_all_networks.md`.**

This file is the how-to.

## Prerequisites

- R with `bnlearn`, `dplyr`, `readr`, `tidyr`.
  - `bnlearn` needs R ≥ 4.4.1 for the current CRAN build. On R 4.4.0 install an
    older source build (e.g. `bnlearn_4.9.4`); if the compile fails with
    `library 'gfortran' not found`, point R at Homebrew's gfortran via
    `~/.R/Makevars`, e.g.
    `FLIBS=-L/opt/homebrew/Cellar/gcc/<ver>/lib/gcc/current -lgfortran -lquadmath`.
- Fitted models under `bn_learning/output/models/` and `bn_cv/output/models/`,
  with their node-mapping CSVs under the matching `output/mappings/`.
- `data/processed/{bn_ready_baseline.csv, flu_response_validation.csv}`.

## Run

From the **repo root** (or from inside `validation/` — scripts locate the repo
root themselves):

```r
Rscript validation/run_all.R          # all five steps in order
# or individually:
Rscript validation/01_assemble_validation_data.R
Rscript validation/02_empirical_benchmark.R
Rscript validation/03_bn_leaf_queries.R
Rscript validation/04_mediation_check.R
Rscript validation/05_validate_all_networks.R   # all networks × {bootstrap, cv}
```

## Layout

```
validation/
  _paths.R                     repo-root finder + shared path anchors (01–04)
  lib_validation.R             shared machinery, source-agnostic (bootstrap|cv)
  01–04_*.R                    L5 deep-dive (bootstrap)
  05_validate_all_networks.R   all networks × both sources
  run_all.R                    orchestrator
  RESULTS.md                   L5 findings
  RESULTS_all_networks.md      multi-network findings
  output/                      all generated CSV/RDS
```

## What each step does

| Script | Purpose | Key output |
|--------|---------|------------|
| `01_assemble_validation_data.R` | Recover subject IDs for the 88 modeled rows (exact feature match — no dependence on the missing CMV-backfill file), join held-out flu targets + antigen summaries | `L5_validation_joined.{rds,csv}` |
| `02_empirical_benchmark.R` | Empirical ground truth: crude + baseline-adjusted root effects, per-antigen breakdown, literature scorecard | `emp_*.csv` |
| `03_bn_leaf_queries.R` | Attach the held-out outcome as a leaf of the roots on the fixed L5 DAG, query `P(seroconvert \| age, cmv)` by likelihood weighting | `bn_leaf_cpqueries.csv`, `bn_leaf_contrasts.csv` |
| `04_mediation_check.R` | Structural test: compose learned `root → cell` edges with held-out `cell → flu`; mediation attenuation | `mediation_edges.csv`, `mediation_attenuation.csv` |
| `lib_validation.R` | Shared machinery: universal subject recovery, flu summaries, leaf query, structural scoring, attenuation — network- **and source-agnostic** (`resolve_files(net, source)`) | (library) |
| `05_validate_all_networks.R` | Recovery + structural validation for **all 9 networks × {bootstrap, cv}**; combined scorecard + flagship axes + attenuation | `all_networks_*.csv`, `leaf_root_benchmark.csv`, `<net>_<source>_validation_joined.rds` |

## Design in one paragraph

The flu response is not a node in any network. Attaching it as a leaf of the
roots (step 3) gives the literal `P(response | roots)`, but that is
**network-invariant** — the internal DAG is irrelevant once the leaf hangs off the
roots. The decisive, per-network test is the **structural** one: each network's
`root → feature` edges were learned without ever seeing flu response, so composing
them with the independently observed held-out `feature → flu` relationship — and
recovering the literature `root → response` direction — is the real evidence the
structure captures biology. Running it under both the bootstrap and the
cross-validated pipelines shows the finding is not an artifact of one
structure-learning method. See `RESULTS_all_networks.md`.
