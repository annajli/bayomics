# Preprocessing Pipeline — Sound Life Multi-Omics → Bayesian Network

**Project:** BPGbio Capstone — Bayesian Network Foundation Model (Sound Life Cohort)
**Cohort:** Allen Institute for Immunology, Sound Life longitudinal healthy-adult cohort
(Gong et al., 2025, *Nature*). Two age groups (Young / Older Adult), sampled across two flu
seasons, with paired multi-omics and clinical measurements.

This directory contains the **Python preprocessing pipeline** that turns raw multi-omics files
into a single, analysis-ready feature matrix for **Bayesian network (BN) structure learning** in
R (`bnlearn`). Unlike neural single-cell foundation models that emit black-box embeddings, the BN
produces an **interpretable joint probability distribution** — a directed acyclic graph over
clinical and molecular features that supports causal-style conditional queries.

---

## What this pipeline produces

Each modality is reduced to a clean, **sample-level, wide-format** block (one row per sample,
biologically selected feature columns). All blocks join on a shared key into one matrix that is
handed to `bnlearn`. A held-out flu-vaccine-response block is produced separately to **benchmark**
the fitted network.

| Deliverable | File | Shape | Role |
|---|---|---|---|
| **BN-ready matrix** | `data/processed/bn_ready_baseline.csv` | 92 × 869 | Training scaffold (all modalities joined) |
| Clinical spine | `data/processed/clinical_baseline_wide.csv` | 92 × 163 | Phenotype layer + discrete roots + canonical metadata |
| Plasma proteomics | `data/processed/olink_baseline_wide.csv` | 92 × 75 | 35 literature-selected proteins |
| Cell-type frequencies | `data/processed/celltype_freq_baseline_wide.csv` | 92 × 123 | Immune composition (CLR) |
| Pseudobulk expression | `data/processed/pseudobulk_baseline_wide.csv` | 92 × 457 | Per-cell-type pathway activity |
| Whole-blood expression | `data/processed/wholeblood_baseline_wide.csv` | 91 × 53 | Hallmark pathway activity |
| Edge constraints | `data/processed/bn_constraints.csv`, `de_whitelist_candidates.csv` | — | Prior-biology whitelist/blacklist for `bnlearn` |
| Feature manifest | `data/processed/feature_manifest.csv` | 1,090 rows | Auditable per-node inclusion rationale |
| **Held-out benchmark** | `data/processed/flu_response_validation.csv` | 92 × 179 | Flu-vaccine response; **never a training input** |

---

## Design principles

These invariants are enforced across every modality so the assembled matrix is coherent and the
prior is defensible to reviewers.

1. **Observations are samples, not cells.** Every modality reduces to one row per
   `sample.sampleKitGuid`. Single-cell data is aggregated to the sample before it enters the BN —
   cells from one sample are pseudoreplicates, not independent observations.
2. **One join key: `sample.sampleKitGuid`**, carried with `subject.subjectGuid`. Preserved in
   every output; never dropped.
3. **Baseline = Flu Year 1, Day 0**, anchored on **`sample.visitName == "Flu Year 1 Day 0"`**
   (92 samples) — *not* `daysSinceFirstVisit == 0`. For ~1/3 of subjects the first-ever draw is an
   earlier visit, so those two definitions disagree; `visitName` is the anchor every block uses.
4. **Informed node selection, structurally flat prior.** Which variables become nodes is chosen
   from **prior biology and literature**, never from this dataset's own variance or association
   structure. We impose *no* prior on which edges exist or their direction — the network learns
   structure freely. Every node's rationale is logged in `feature_manifest.csv`, so the prior is
   auditable and demonstrably non-circular.
5. **No imputation during preprocessing.** Missingness is measured and flagged (features > 20%
   missing are marked `exclude_candidate`); imputation is a downstream group decision.
6. **Raw data is immutable.** Nothing in `data/raw/` is overwritten. All outputs go to
   `data/processed/`.
7. **Flag, don't drop or overwrite.** Quality issues ride alongside each value as paired boolean
   `*_flag` columns; only truly empty (100%-missing) columns are removed.

---

## Pipeline flow

```
        ┌──────────────────────────────────────────────┐
        │  STEP 1 — CLINICAL SPINE                       │
        │  Clinical labs + canonical subject/sample      │
        │  metadata, filtered to baseline (92 samples)   │
        │  → defines the sample set everyone joins to    │
        └───────────────────────┬────────────────────────┘
                                │  (join on sample.sampleKitGuid)
   ┌──────────────┬─────────────┼──────────────┬───────────────┐
   ▼              ▼             ▼              ▼               ▼
 Olink        Cell-type     Pseudobulk    Whole-blood     Differential
 proteomics   frequencies   expression    RNA-seq         expression
                                                          (→ edge constraints)
        STEP 2 — MOLECULAR BLOCKS (independent, parallel)
                                │
                                ▼
        ┌──────────────────────────────────────────────┐
        │  STEP 3 — JOIN → bn_ready_baseline.csv         │
        └───────────────────────┬────────────────────────┘
                                ▼
        ┌──────────────────────────────────────────────┐
        │  VALIDATION (produced independently, held out) │
        │  Flu serology & HAI → flu_response_validation  │
        └──────────────────────────────────────────────┘
```

The **clinical spine must exist first** (every molecular block joins against it). The molecular
blocks are mutually independent. Flu-vaccine response is a **held-out benchmark**, produced
independently and never fed into training.

---

## Modality summaries

### 1. Clinical labs & metadata (the spine) — `clinical_labs_metadata_preprocessing.ipynb`
Source: `sound_life_labs_metadata.csv` + `sound_life_clinical_lab_descriptors.csv`.
The source of truth for `subject.*` / `sample.*` metadata and the baseline sample set.

- **53 clinical phenotype nodes** (blood chemistry, CBC/differential, lipids, anthropometrics,
  inflammatory markers) + derived **NLR** (neutrophil-to-lymphocyte ratio, an inflammaging marker).
- **3 discrete root nodes:** age group, biological sex, CMV serostatus — the upstream drivers of
  immune variation. Modeled as roots because they sit upstream of essentially all immune biology.
- **Two independent QC flag families**, kept deliberately separate: `_flag_outOfRefRange`
  (outside the clinical *normal* range — clinically abnormal but **not** a data error) vs
  `_flag_implausible` (physically impossible values — genuine data errors). Conflating them would
  mislabel real age/CMV biology as bad data.
- The clinical panel is an already-curated standard-of-care assay set, so it is kept whole rather
  than variance-filtered.

### 2. Plasma proteomics (Olink) — `olink_preprocessing.ipynb`
Source: `sound_life_all_olink.csv` (~1,470 proteins across 4 Olink Explore panels).

- Reduced to a **35-protein panel (14 core + 21 extended)**, each with an **independent,
  non–Sound-Life immune-aging citation** — protecting the flat prior (selection is external
  biology, never this dataset's variance). Full evidence log in `context/olink_panel_selection.md`.
- Value column `NPX_norm` (batch-normalized). Below–limit-of-detection values flagged per protein
  (`_flag_belowLOD`); a candidate is disqualified if below LOD in > 30% of baseline samples.
- Bridging proteins (IL6/TNF/CXCL8, measured on all panels) resolved by fixed panel priority —
  never averaged, since panels carry systematic offsets.

### 3. Cell-type frequencies — `cell_type_freq_preprocessing.ipynb`
Source: `sound_life_AIFI_L{1,2,3}_frequencies.csv` (already sample-level).

- All three annotation levels retained (L1 = 9, L2 = 29, L3 = 83 = **121 composition features**);
  the recommended BN subset is **L2 (29 intermediate cell types)** — interpretable and within
  budget. The finer/coarser levels are carried for flexibility and subset at fit time.
- Value = **pseudocount-adjusted CLR** (`AIFI_Lx_clr_pseudo`), the correct continuous
  representation for compositional data.

### 4. Pseudobulk scRNA-seq — `pseudobulk_preprocessing.ipynb`
Source: `sound-life_AIFI_L3_pseudobulk.h5ad` (cell-type × sample).

- Collapsed to **per-cell-type pathway activity**: **65 retained L3 cell types** × **6 immune
  Hallmark pathways** (interferon-α/γ, inflammatory, TNFα-NFκB, IL2-STAT5, IL6-JAK-STAT3) =
  **390 module scores**, scored as **mean z-score within cell type** across samples. (Cell types
  too sparse across the cohort are dropped rather than imputed.)
- Per-cell-type reliability flag (`_low_ncells_flag`, 65 columns): scores set missing when a cell
  type has too few cells in a sample. Genes are never used individually — pathway modules keep
  dimensionality tractable at n = 92.

### 5. Whole-blood RNA-seq — `whole_blood_preprocessing.ipynb`
Source: `sound-life_whole-blood_no-stim.h5ad` (already one sample per row).

- CPM-normalized, log-transformed, then reduced to **50 MSigDB Hallmark pathway scores via
  ssGSEA** (within-sample ranking, appropriate for already-sample-level data).
- Library-size QC flag (`wb.low_lib_flag`); low-library samples have scores set missing.
- **91 samples** (one baseline sample absent from this assay) — the expected cross-modality
  coverage gap, reconciled at join time rather than imputed.

### 6. Differential-expression constraints — `differential_expression_preprocessing.ipynb`
Source: the DESeq2 contrast tables (age group, sex, CMV, flu vaccine).

- These are **summary statistics, not observations** — never rows or columns in the BN matrix.
- Translated into **edge constraints** for `bnlearn`: `bn_constraints.csv` (blacklist — e.g.
  molecular nodes may not point *into* the discrete roots) and `de_whitelist_candidates.csv`
  (biology-first whitelist, e.g. CMV → interferon-γ pathways, each cited). Constrains edge
  *direction* from prior biology without selecting nodes from the same data.

### 7. Join & assemble — `bn_ready_baseline_join.ipynb`
- **Left-joins** each molecular block onto the clinical spine on `sample.sampleKitGuid`. The spine
  is a strict superset (no orphan samples), guaranteeing the **92-row** invariant.
- **Full assembly** (~652 candidate node columns): node *selection* stays downstream in R, so the
  flat prior is preserved and each modality owner retains selection authority.
- `coverage.*` columns record per-sample block presence (**91 of 92 samples have all modalities**;
  the 1 partial lacks whole-blood, left missing — no imputation). All `*_flag` columns are ordered
  last.

### 8. Flu serology & HAI — held-out validation — `flu_validation_preprocessing.ipynb`
Source: `sound-life_flu_serology_single.csv`, `sound-life_flu_hai_single.csv`.

- **Not a training input.** Serology is measured only at vaccine visits, so including it as a node
  would inject structured missingness. It exists purely to **benchmark** the fitted network.
- Per-antigen response targets: **IgG fold change** `conc(Day 7)/conc(Day 0)`, **seroconversion**
  (≥ 4-fold rise), and **peak HAI inhibition** (Day 7). Pre-existing-immunity flags mark high
  Day-0 titers (ceiling effects that suppress the Day-7 response).
- One row per subject, keyed on the Flu-Year-1 Day-0 sample so it aligns to the fitted-network
  subjects at benchmark time. IgG and HAI panels kept as separate blocks (their antigen panels
  differ). HAI plate replicates are averaged, with a disagreement flag where plates diverge.
- **How it's used:** once the network is fit, query it for conditional probabilities against
  literature-established expectations — e.g. *P(seroconversion | Older Adult, CMV⁺)*, where age and
  CMV are known to suppress vaccine response. Agreement between the network's inference and the
  established direction/magnitude is evidence the learned structure captures real biology.

---

## Repository layout

```
preprocessing/
  clinical_labs_metadata_preprocessing.ipynb   # spine (Step 1)
  olink_preprocessing.ipynb                     # molecular blocks (Step 2)
  cell_type_freq_preprocessing.ipynb
  pseudobulk_preprocessing.ipynb
  whole_blood_preprocessing.ipynb
  differential_expression_preprocessing.ipynb   # edge constraints
  bn_ready_baseline_join.ipynb                  # join (Step 3)
  flu_validation_preprocessing.ipynb            # held-out benchmark
src/
  load.py          # GUID-safe CSV loading
  qc.py            # duplicate / missingness checks
  transforms.py    # shared transforms (snake-case, CLR/score helpers)
  io.py            # export-contract enforcement, manifest append
data/
  raw/             # immutable source files
  processed/       # all pipeline outputs (see table above)
context/           # design docs, panel-selection evidence logs
```

Each preprocessing notebook follows the same shape: **ingest → baseline filter → QC / missingness
audit → modality-specific transform → biological feature selection → export**, and writes both its
wide block and its rows in `feature_manifest.csv`.

---

## Reproducibility

- **Environment:** `ds6001` conda environment (Python 3.11+; `pandas`, `numpy`, `scipy`, `scanpy`,
  `anndata`, `gseapy`, `scikit-learn`, `missingno`, `matplotlib`, `seaborn`).
- **Run a notebook:**
  ```bash
  conda run -n ds6001 jupyter nbconvert --to notebook --execute --inplace \
    preprocessing/<notebook>.ipynb
  ```
- **Order:** the clinical spine must be built before the join; molecular blocks and the flu
  validation block can run in any order. Every notebook is self-checking (asserts keys, row counts,
  and the export contract) and re-runnable.

---

## Handoff to the modeling phase

`bn_ready_baseline.csv` and `bn_constraints.csv` feed the **R / `bnlearn`** phase: structure
learning and parameter fitting under a **conditional linear Gaussian (CLG)** model (continuous
clinical/omics nodes with discrete age-group / sex / CMV roots), seeded by the prior-biology edge
constraints. The fitted network is then benchmarked against `flu_response_validation.csv`.

## Notes for the sponsor

- **Small n, wide p by design.** Baseline n ≈ 92. The assembled matrix is intentionally wide;
  final node selection (targeting ~80–120 nodes) happens at fit time in R, keeping every selection
  decision biologically motivated and auditable.
- **The prior is defensible.** Node inclusion is driven by external literature (logged per node in
  `feature_manifest.csv`), and edge constraints come from independent differential-expression
  biology — the network is not told which edges to find. This is what makes the resulting DAG an
  interpretable, non-circular model rather than a data-dredged one.
- **Validation is genuinely held out.** Flu-vaccine response never touches training, so
  benchmarking against it is an honest test of whether the learned structure recovers known
  immunology (age and CMV effects on vaccine response).
