# bn_cv Pipeline — Complete Results Summary

## What each network is

| Network | Modality | Description | Node source |
|---|---|---|---|
| **L1a** | Clinical | Full clinical panel — every clinical lab/measurement column not excluded for missingness or determinism | `config/L1a_clinical_full.R` |
| **L1b** | Clinical | Curated clinical panel — a literature-selected subset of ~22 clinically meaningful markers, smaller and more sample-efficient than L1a | `config/L1b_clinical_curated.R` |
| **L2** | Olink | Proteomics — 34-protein inflammation/immune panel measured via Olink NPX assay | `config/L2_olink_adapted.R` |
| **L3a** | Combined | L1a (full clinical) + L2 (Olink) merged into one joint network | `config/L3a_clinical_full_olink.R` |
| **L3b** | Combined | L1b (curated clinical) + L2 (Olink) merged into one joint network | `config/L3b_clinical_curated_olink.R` |
| **L4** | Whole blood | Bulk RNA-seq scored against the 50 MSigDB Hallmark gene sets — pathway-level activity summaries, not individual gene expression | `config/L4_wb_pathways.R` |
| **L5** | Immunophenotyping | Flow/mass cytometry cell-type frequencies, in a 3-level gating hierarchy (L1: 9 broad lineages, L2: 29 subtypes, L3: 71 finest subtypes after dedup) | `config/L5_freq_immunophenotype.R` |
| **L6** | PBMC signaling | Phospho-flow cytometry — PBMCs stimulated with 6 signals, measured for intracellular pathway activation across 24 cell types (144 = 24 × 6 nodes) | `config/L6_pb_signaling.R` |

All networks share the same 3 root nodes: `age_group`, `sex`, `cmv` (CMV serostatus).

Every network except L5 uses the shared `01_data_prep_cv.Rmd` for data
preparation. **L5 uses a dedicated `01_data_prep_L5_cv.Rmd`**, since it's
the only network needing the sibling log-ratio (ALR) transform described
below. All networks use the shared `02_structure_learning_cv.Rmd` for
structure learning.

---

## Cross-cutting fixes (apply across some or all networks)

| Fix | What it does | Applies to |
|---|---|---|
| Empirical-valley conservative threshold | Replaced bootstrap's fixed 0.85 cutoff with a per-network histogram-based valley in the CV edge-strength distribution — CV strengths cluster near 0/1 much more than bootstrap's, so a literal 0.85 doesn't mean "conservative" the same way | All networks |
| WBC-cluster blacklist (Rule 3) | Forbids direct edges between WBC and its differential white-cell-count components (defensively scoped via `intersect()` — no-ops where the cluster isn't present) | L1a (full 6-node cluster), L1b (5 nodes, no basophils); no-ops elsewhere |
| Hierarchy blacklist (Rule 4) | Forbids direct parent-child edges using a config-supplied `PARENT_MAP` | L5 only; no-ops for all others |
| `maxp = 10` | Caps parent count per node in `hc()`/`tabu()`/`mmhc()`/`bn.cv()`'s internal search — prevents nodes from accumulating enough parents to produce near-singular local fits | All networks (confirmed present and re-verified via rerun in all 6 working networks) |

---

## L1a — Clinical (full panel)

**37 nodes** (3 roots + 34 continuous), **81/92 samples**.

**Fixes applied:**
- 10 deterministic/derived columns excluded (MCH/MCHC/MCV, 5 `perc_*`
  columns, `non_hdl`, `hdl_ratio`) — each is an exact arithmetic function of
  other columns already in the network (e.g. MCHC = Hgb/Hct×100), the same
  category of problem as the already-excluded NLR.
- `ldl`/`total_chol` checked against the same concern and explicitly kept —
  verified the relationship isn't a closed formula (residual spread up to
  23.4 mg/dL, correlates with triglyceride level, consistent with a
  variable-divisor lab calculation rather than the fixed Friedewald formula
  these 3 columns alone would need to reproduce it exactly).
- WBC-cluster blacklist rule added after this network repeatedly threw
  "would introduce cycles" warnings — WBC ≈ sum of its 5 differential
  counts, all independently present as nodes.

**Results (maxp-confirmed via rerun):**
```
Tabu: 80 edges, BIC-cg = -6485.19
MMHC: 29 edges, BIC-cg = -6676.35
avg_opt (t=0.5): 67 edges
avg_conservative (t=0.625): 55 edges
Held-out loss: mean = 85.115, sd = 5.679  (sd/mean ≈ 7%)
```
`maxp` changed nothing here (ratio 37/81 ≈ 0.46, well below where the
constraint binds) — numbers are bit-identical pre/post.

**History worth knowing:** this network's *original* CV run (before the
deterministic-column fix) had a catastrophic loss (mean in the hundreds of
millions, sd ~10x the mean, traced to a single dominant fold). Removing the
10 deterministic columns resolved nearly all of that instability on its
own — suggesting the original "too many nodes for the sample size"
hypothesis was likely secondary to collinearity for this network.

---

## L1b — Clinical (curated panel)

**25 nodes** (3 roots + 22 continuous), **87/92 samples**.

**Fixes applied:** WBC-cluster blacklist (partial — this panel excludes
`basophils`, so the rule correctly scopes to the 5 remaining members).
**Open item:** this panel was never reviewed for the broader
deterministic-column issue found in L1a — only NLR is excluded here.

**Results (maxp-confirmed via rerun):**
```
Tabu: 49 edges, BIC-cg = -4955.63
MMHC: 15 edges, BIC-cg = -5044.08
avg_opt (t=0.5): 40 edges
avg_conservative (t=0.675): 29 edges
Held-out loss: mean = 58.867, sd = 2.781  (tightest CV loss of all 6 working networks)
```
`maxp` changed nothing here (ratio 25/87 ≈ 0.29, the lowest of any network).

---

## L2 — Olink proteomics

**37 nodes** (3 roots + 34 continuous), **88/92 samples**.

**Fixes applied:** none needed beyond the pre-existing MMP9 missingness
exclusion. No deterministic-collinearity issues found or expected — Olink
NPX values don't have clinical data's formula-derived-column problem.

**Results (maxp-confirmed via rerun):**
```
Tabu: 91 edges, BIC-cg = -2368.16
MMHC: 36 edges, BIC-cg = -2526.07
avg_opt (t=0.5): 66 edges
avg_conservative (t=0.575): 60 edges
Held-out loss: mean = 32.390, sd = 6.376  (sd/mean ≈ 20%, flattest/least
                                            sharply-defined threshold valley
                                            of the clinical-scale networks,
                                            not investigated further —
                                            still a normal loss range)
```
`maxp` changed nothing here (ratio 37/88 ≈ 0.42).

---

## L3a — Combined: Clinical (full) + Olink

**71 nodes** (3 roots + 34 clinical + 34 Olink), config built by merging
L1a's and L2's already-fixed node sets and exclusion lists (not re-derived
independently).

**STATUS: BROKEN, UNRESOLVED.** Rendered once — showed a catastrophic CV
loss (same broad symptom pattern as L5's original problem: high sd/mean
ratio). **Never diagnosed or fixed** — work moved to other networks before
this was resolved. Do not treat L3a's edges/thresholds as trustworthy.
Needs the same diagnostic process used for L5 (sort fold losses → identify
dominant fold(s) → check for near-singular node fits → likely `maxp` or a
cross-modality collinearity issue not yet checked).

**Known gap, separate from the above:** neither L1a's nor L2's individual
review checked for collinearity *between* the two panels (e.g. a clinical
inflammatory marker highly correlated with a specific Olink cytokine for a
real biological reason) — this is a distinct, unexplored possible
contributor.

---

## L3b — Combined: Clinical (curated) + Olink

**59 nodes** (3 roots + 22 clinical + 34 Olink), config built the same way
as L3a but from L1b's curated panel.

**STATUS: NEVER RENDERED.** Config exists and was verified for
alias-collision safety at creation time, but this network has not been run
through either `01_data_prep_cv.Rmd` or `02_structure_learning_cv.Rmd` at
all in this project.

---

## L4 — Whole blood (Hallmark pathway scores)

**53 nodes** (3 roots + 50 continuous), **86/92 samples**.

**Fixes applied:**
- `wb.low_lib_flag` excluded from `CONTINUOUS_MAP` — a boolean RNA-seq QC
  flag (90 False / 1 True / 1 NA), not a pathway score. Caught before it
  was ever included as a continuous node.
- No deterministic-collinearity issues found among the 50 real pathway
  scores (Hallmark sets can share member genes and correlate biologically,
  but aren't exact arithmetic functions of one another).

**Results (maxp-confirmed via rerun):**
```
Tabu: 139 edges, BIC-cg = 11738.64   <- positive BIC-cg, see note below
MMHC: 46 edges, BIC-cg = 11338.03
avg_opt (t=0.5): 97 edges
avg_conservative (t=0.575): 89 edges
Held-out loss: mean = -131.685, sd = 3.632  (tightest CV loss, in
                                              absolute-value terms, of any
                                              network in the project)
```
`maxp` changed nothing here (ratio 53/86 ≈ 0.62).

**On the positive BIC-cg:** this same pattern recurred in L6 and was
directly investigated there (per-node residual sd and parent-count check —
see L6 section). Same explanation likely applies here (BIC's complexity
penalty behaving atypically on the single full-data fit at this
node-to-sample ratio, not hidden numerical instability), though the L6-style
node-level check was not separately repeated for L4 — its CV loss was
already excellent, so this was lower priority.

---

## L5 — Immune cell frequencies (3-level gating hierarchy)

**91 nodes** (3 roots + 88 continuous, after the sibling log-ratio
transform), **88/92 samples**. The most heavily-debugged network in the
project — original data was 121 raw columns across 3 nested levels
(L1: 9 broad lineages, L2: 29 subtypes, L3: 83 finest subtypes).

**Fixes applied, in the order they were found:**

1. **12 exact-duplicate columns removed.** 6 base column names each had 4
   suffixed variants (`base`, `.1`, `.2`, `.3`) — confirmed `base == .2`
   and `.1 == .3` exactly, in all 92 subjects, across all 6 groups, no
   exceptions (almost certainly a source-data concatenation artifact). Kept
   `base` and `.1` (renamed `_alt`), dropped `.2`/`.3`.
2. **Hierarchy blacklist (Rule 4)** — forbids direct edges between a node
   and its immediate gating-tree parent, via a `PARENT_MAP` built into the
   config. Needed a defensive `intersect()` guard for networks where the
   cluster isn't fully present.
3. **Two real implementation bugs**, same failure pattern both times: a
   blacklist component (`bl_hierarchy`, then separately `bl_wbc`) was
   computed and its summary correctly printed, but never actually added
   into the final `blacklist <- rbind(...)` call — silently not applied
   despite looking like it was. Caught by cross-checking the printed row
   count against the actual total blacklist size.
4. **Additive log-ratio (ALR) transform.** Direct parent-child blacklisting
   alone wasn't sufficient — siblings under the same parent are also
   compositionally redundant (verified: correlation between
   `Σexp(children)` and `exp(parent)` ≈ 0.87, real but partial redundancy,
   not an exact formula). For each sibling group, replaced all but the
   highest-mean-abundance sibling with `log(sibling) - log(reference)`, and
   dropped the reference node. Extended beyond `PARENT_MAP`'s explicit
   scope to also treat the 9 L1 categories as one additional top-level
   sibling group (mutually compositional — should sum to ~100% of
   measured cells — but had no parent entry of their own).
5. **`maxp = 10`** — the actual root cause of the remaining catastrophic
   instability after the above fixes reduced but didn't eliminate it.
   Diagnosed by: sorting fold losses (top 3 of 100 folds accounted for
   >90% of total loss) → finding the worst fold's fitted model had a node
   (`l1_progenitor_cell_alr`) with residual sd = 3.13e-07 (near-singular) →
   tracing this to that node having 31 parents in an unconstrained search
   with only ~79 training samples per fold. Capping `maxp` fixed it
   directly.

**Results (final, post all fixes, re-verified):**
```
Tabu: 431 edges, BIC-cg = -3637.73
MMHC: 90 edges, BIC-cg = -6743.73
avg_opt (t=0.5): 297 edges
avg_conservative (t=0.875): 130 edges
Held-out loss: mean = 61.211, sd = 15.910  (sd/mean ≈ 26% — loss histogram
                                             shows a real right tail out to
                                             ~120, consistent with this
                                             ratio, not a hidden problem)
```

Progression of the CV loss mean across each fix stage, for reference:
7.4 billion (blacklist bug) → 715 million (blacklist fixed) → 570 million
(ALR added) → **61.2 (maxp added)** — ~8 orders of magnitude improvement,
only the last step of which was sufficient on its own; the earlier steps
were real, necessary fixes too.

**Post-`maxp` verification:** directly checked `fit_tabu`'s per-node
residual sd and parent counts — no near-singular fits remain (smallest sd
values are all ordinary, nowhere near the earlier 3.13e-07 failure).

---

## L6 — PBMC functional/signaling response

**147 nodes** (3 roots + 144 continuous), **77/92 samples** — the highest
node-to-sample ratio of any network in this project (≈1.9).

**Fixes applied:**
1. **Cell-type selection restricted to 24 of 65 measured types.** The
   other 41 were dropped entirely, based on missingness driven by
   cell-type rarity (if a subject didn't have enough cells of a given
   type, all 6 signaling readouts for that type go missing together —
   confirmed 100% consistent within-type missingness pattern). A 5%
   per-cell-type missingness cutoff was chosen to stay under the project's
   existing 15% total-sample-loss precedent (used earlier to justify the
   bmi/height/MMP9 exclusions) — keeps complete-case n at 77/92 vs.
   crashing to 18/92 at a 20% cutoff.
2. **65 `*_low_ncells_flag` QC columns excluded** — same category as
   `wb.low_lib_flag`, boolean metadata not biology.
3. **`maxp = 10`** applied from the start, given the lesson already
   learned from L5.

**Checked, no exclusion needed:** pairwise correlations between the 6
stimulation conditions within each of the 24 kept cell types (360 pairs).
46 pairs exceeded r=0.9, overwhelmingly `ifn_alpha` vs. `ifn_gamma` (~20/24
cell types) — consistent with both interferon types converging on shared
STAT1 phosphorylation signaling, real biology, not an arithmetic identity
(max correlation 0.976, well short of the near-1.0 exact-duplicate
threshold that justified exclusions elsewhere). Left unblacklisted, same
standard applied to the clinical electrolyte cluster earlier in the
project.

**Results:**
```
Tabu: 721 edges, BIC-cg = 10238.35   <- positive BIC-cg
MMHC: 141 edges, BIC-cg = 6671.32
avg_opt (t=0.5): 349 edges
avg_conservative (t=0.875): 145 edges
Held-out loss: mean = -83.477, sd = 15.445  (broad but smooth distribution,
                                              no outlier spike)
```

**On the positive BIC-cg — directly investigated, not assumed benign:**
pulled `fit_tabu`'s per-node residual sd and parent counts directly.
Smallest sd values were all in the 0.04-0.05 range (~100,000x larger than
L5's original 3.13e-07 failure signature). Several nodes sit exactly at the
`maxp=10` ceiling, confirming the cap is actively constraining the search.
The lowest-sd node was also at the parent-count ceiling, and its sd (0.040)
was still entirely ordinary. Conclusion: the positive BIC-cg reflects BIC's
complexity penalty behaving atypically at this node-to-sample ratio on the
single full-data fit, not a hidden numerical problem.

---

## Summary table

| Network | Nodes | Samples | Ratio | Tabu edges | avg_conservative | CV loss mean | CV loss sd | Status |
|---|---|---|---|---|---|---|---|---|
| L1a | 37 | 81 | 0.46 | 80 | 55 (t=0.625) | 85.1 | 5.7 | ✅ Verified |
| L1b | 25 | 87 | 0.29 | 49 | 29 (t=0.675) | 58.9 | 2.8 | ✅ Verified |
| L2 | 37 | 88 | 0.42 | 91 | 60 (t=0.575) | 32.4 | 6.4 | ✅ Verified |
| L3a | 71 | ~75 | ~1.05 | — | — | catastrophic | catastrophic | ❌ Broken, unresolved |
| L3b | 59 | ~79 | ~0.75 | — | — | — | — | ⚠️ Never rendered |
| L4 | 53 | 86 | 0.62 | 139 | 89 (t=0.575) | -131.7 | 3.6 | ✅ Verified |
| L5 | 91 | 88 | 1.03 | 431 | 130 (t=0.875) | 61.2 | 15.9 | ✅ Verified |
| L6 | 147 | 77 | 1.9 | 721 | 145 (t=0.875) | -83.5 | 15.4 | ✅ Verified |

All six working networks show healthy, non-catastrophic CV loss, confirmed
both numerically (console output) and visually (edge-strength and loss
histograms all show the expected shape — signal/noise separation in edge
strength, smooth single-cluster loss distributions with no outlier spike).
`maxp=10` is confirmed present and re-verified via full rerun in all six;
it changed results only for L5 and L6 (the two networks with ratio > 1.0),
exactly as predicted from the node-to-sample ratio reasoning.

