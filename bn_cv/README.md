# bn_adapted CV Pipeline — Session Changelog

This documents all fixes and changes made to the `bn_adapted` (cross-validation)
pipeline in this session, why each was needed, and what's still open. Scope is
limited to `bn_adapted/` — the original `bn_learning/` (bootstrap) pipeline was
**not** touched by anything below.

---

## Files changed

| File | What changed |
|---|---|
| `02_structure_learning_cv.Rmd` | Conservative threshold logic rewritten; WBC-cluster blacklist rule added |
| `03_model_comparison_cv.Rmd` | Multiple filename-pattern bugs fixed; `avg_85` references updated; Section 5 restored |
| `config/L1a_clinical_full_v2.R` | 10 deterministic/derived columns removed from `CONTINUOUS_MAP` |

`L1b_clinical_curated_v2.R` and `L2_olink_v2.R` were **not edited** — see
"Not changed" section below.

---

## 1. `02_structure_learning_cv.Rmd`

### 1a. Conservative threshold — replaced fixed 0.85 with an empirical valley

**Problem:** the original conservative threshold was hardcoded to 0.85,
copied from the bootstrap pipeline's convention. CV-derived edge strengths
cluster much more heavily near 0 and 1 than bootstrap strengths do (CV folds
share ~90% of their training data, so an edge either survives nearly every
fold or almost none), so a threshold value that means "conservative" under
bootstrap doesn't mean the same thing under CV. Applying 0.85 literally gave
inconsistent, not-really-comparable edge counts across networks (e.g. L1b:
38→10 edges under bootstrap's 0.85 vs. 42→24 edges under a literal 0.85 on
CV strengths).

**Fix, in two iterations:**
1. First attempt used `density()` (KDE) to find the local minimum ("valley")
   in the mid-range (0.3–0.95) of the strength distribution. This worked for
   L1b but failed for L1a — the KDE got pulled toward the search boundary by
   L1a's large near-1.0 edge cluster, returning t=0.95 (the boundary itself)
   instead of the true valley at t≈0.82.
2. **Final fix:** replaced KDE with simple histogram binning (0.05-wide bins
   over 0.3–0.95), taking the bin with the minimum count as the threshold.
   Binning can't overshoot a boundary the way a smoothed density estimate can.
   Also prints the raw bin counts to the render log so the valley can be
   eyeballed against the edge-strength histogram plot.

**Variable renamed:** `avg_85` → `avg_conservative` throughout (the value is
no longer literally 0.85, so the old name became misleading). This rename is
also reflected in the `models <- list(...)` save block, meaning **every
`*_models_cv.rds` file now stores this network under the key
`avg_conservative`, not `avg_85`** — anything downstream that reads these
files needs to use the new name (see Section 3 below, this is exactly what
broke `03_model_comparison_cv.Rmd`).

**Resulting thresholds per network (this session's final runs):**

| Network | Empirical valley threshold | avg_conservative edges | avg_opt (t=0.5) edges |
|---|---|---|---|
| L1a_clinical_full (post config fix) | 0.625 | 55 | 67 |
| L1b_clinical_curated | 0.675 | 29 | 40 |
| L2_olink | 0.575 | 60 | 66 |

Note L2's valley is shallower/less sharply defined than L1a/L1b's — worth
treating as directionally correct but lower-confidence than the other two.

### 1b. WBC-cluster blacklist rule (Rule 3)

**Problem:** `averaged.network()` threw repeated "would introduce cycles,
ignoring" warnings at loose thresholds (t=0.5) for L1a. Root cause: WBC is
(approximately) the sum of its five differential white cell counts
(neutrophils, lymphocytes, monocytes, eosinophils, basophils), which are
*also* individually in the network — the same deterministic-collinearity
problem that had already gotten NLR excluded, just not caught for WBC.

**Fix:** added a blacklist rule forbidding edges between `wbc` and its
present differential-count components, in both directions:

```r
wbc_components <- c("neutrophils", "lymphocytes", "monocytes", "eosinophils", "basophils")
wbc_cluster <- intersect(c("wbc", wbc_components), all_nodes)

if (length(wbc_cluster) >= 2) {
  bl_wbc <- expand.grid(from = wbc_cluster, to = wbc_cluster, stringsAsFactors = FALSE) %>%
    filter(from != to)
} else {
  bl_wbc <- data.frame(from = character(0), to = character(0))
}

blacklist <- rbind(bl_into_roots, bl_between_roots, bl_wbc) %>% distinct()
```

**Important implementation detail:** the component list is intersected with
`all_nodes` (the actual node set for whichever network is currently loaded
via `params$config`) before building the blacklist. This is required because
not every network config includes `wbc` or all 5 components — L1b's curated
panel omits `basophils`, and L2 (Olink) has none of these nodes at all.
Referencing a blacklist node that doesn't exist in `bn_df` throws a hard
`build.blacklist()` error rather than failing silently, so this defensive
`intersect()` is load-bearing, not just tidy — an earlier hardcoded version
of this rule broke L1b's render with exactly that error.

This chunk is shared across all three networks via `params$config`, so this
fix applies automatically to all of them.

---

## 2. `config/L1a_clinical_full_v2.R`

**Problem:** beyond WBC, L1a's clinical panel contains several other
columns that are directly calculated from other columns already in the
network — the same category of issue as the already-excluded NLR, just not
previously caught. These produced additional cycle warnings
(`mchc -> hemoglobin`, `perc_neutrophils -> neutrophils`, etc.) and,
more importantly, meant the network could contain "edges" that are really
just restated arithmetic rather than discovered relationships — a problem
for a project whose stated goal is an *interpretable* network.

**Columns removed from `CONTINUOUS_MAP`** (44 → 34 entries):

| Column(s) | Why |
|---|---|
| `mch`, `mchc`, `mcv` | Standard CBC indices computed directly from `hemoglobin`, `hematocrit`, `rbc` by the lab analyzer (MCH = Hgb/RBC×10, MCHC = Hgb/Hct×100, MCV = Hct/RBC×10) |
| `perc_basophils`, `perc_eosinophils`, `perc_lymphocytes`, `perc_monocytes`, `perc_neutrophils` | Each is `count / wbc × 100` |
| `non_hdl`, `hdl_ratio` | `non_hdl = total_chol - hdl`; `hdl_ratio = total_chol / hdl` |

**Checked but NOT excluded — `ldl` and `total_chol` both retained:**
verified against `bn_ready_baseline.csv` (n=91 complete cases) that these
are *not* a simple Friedewald calculation of each other. Residual
`total_chol - (hdl + ldl + triglycerides/5)` has real spread (up to 23.4
mg/dL off, correlation 0.993 but not ≈1.0), and the residual grows with
triglyceride level even below the classic Friedewald breakdown point
(TG>400) — consistent with a Martin-Hopkins (variable-divisor) LDL
calculation or a direct assay, neither of which is reproducible from the
three columns in this dataset alone. This is a genuine, if partially
correlated, degree of freedom — unlike NLR/MCHC/percentages, which are
exact closed-form identities.

**Effect on structure learning (L1a, before vs. after this config change):**

| | Before | After |
|---|---|---|
| Nodes | 47 (3 roots + 44 continuous) | 37 (3 roots + 34 continuous) |
| Tabu edges | 132 | 80 |
| MMHC edges | 40 | 29 |
| CV loss mean / sd | 565,915.79 / 5,655,016.09 | **85.115 / 5.679** |

The CV loss collapse (sd was ~10x the mean before, ~7% of the mean after) is
the headline result — strongly suggests the originally-reported "one
catastrophic CV fold" instability (documented in the original project
handoff, attributed there to a high node-to-sample ratio) was actually
being driven substantially by these deterministic/collinear columns rather
than sample size alone. Worth revisiting that original framing.

**Documentation-only, no functional effect:** `DETERMINISTIC_EXCLUDE` was
also updated to list all 10 newly-excluded columns for documentation
purposes, matching the existing NLR entry. Per the project's own prior
finding, **this list is never read by any code** — only a column's absence
from `CONTINUOUS_MAP` actually excludes it. This was caught during this
session (an early draft of this fix only updated `DETERMINISTIC_EXCLUDE`
and would have been a silent no-op) — flagging again here since it's a
recurring failure mode on this project.

---

## 3. `03_model_comparison_cv.Rmd`

This file needed the most fixes, all in the same family: filenames that had
drifted to include a `_cv` suffix (per your naming decision to keep CV
output from colliding with bootstrap output), and the `avg_85` →
`avg_conservative` rename from Section 1 above, without every reference in
this file being updated to match. Fixed in the order they were hit:

| Location | Was | Now |
|---|---|---|
| Section 1 (`load-models`) | `pattern = "_models\\.rds$"` | `pattern = "_models_cv\\.rds$"` (and matching `gsub()`) |
| Section 2 (`summary-table`) | `narcs(m$avg_85)` | `narcs(m$avg_conservative)` |
| Section 2 (`summary-table`) | `readRDS(sprintf("output/models/%s_data.rds", nm))` | `readRDS(sprintf("output/models/%s_data_cv.rds", nm))` |
| Section 4 (`sanity-queries`) | `readRDS("output/models/L1b_clinical_curated_data.rds")` | `readRDS("output/models/L1b_clinical_curated_data_cv.rds")` |
| Section 6 (`handoff`) | `pattern = "_node_mapping\\.csv$"` | `pattern = "_node_mapping_cv\\.csv$"` |
| Section 6 (`handoff`) | example `cat()` string referencing `L1b_clinical_curated_models.rds` | updated to `..._models_cv.rds` |

**Section 5 (CV loss boxplot across networks)** was temporarily dropped from
one version of the file mid-session and restored to its correct position
between Sections 4 and 6, content unchanged (it has no filename
dependencies, so it needed no fixes).

**Not changed, checked and confirmed fine:**
- Section 3's hardcoded network names (`L1a_clinical_full`,
  `L1b_clinical_curated`) — these match the names produced by Section 1's
  `gsub()` correctly, since `NETWORK_ID` values were never actually renamed
  with a `_cv` suffix, only the *saved file* names were.
- Section 6's generic `list.files("output/models", pattern = "\\.rds$")` —
  intentionally broad, matches both `_data_cv.rds` and `_models_cv.rds`
  files correctly as-is.

**Final render result (all three networks loaded successfully):**
- Cross-level edge comparison: L1a (80 edges) vs. L1b (49 edges) → 22 shared,
  58 L1a-only, 27 L1b-only.
- All 4 sanity queries ran; 3 of 4 showed the expected direction with a
  clear margin (CRP by age, hemoglobin by sex, CRP by age+CMV jointly — the
  last of these on a narrow ~0.02 margin, worth treating as weak evidence
  rather than a clean pass given n~90). Lymphocyte/CMV query has no
  hardcoded PASS/FAIL check in the original file, just directional output.

---

## Not changed / explicitly out of scope this session

- **`config/L1b_clinical_curated_v2.R`** — not edited. L1b already excludes
  `basophils` from its differential-count panel (curation predates this
  session), and the WBC blacklist rule's `intersect()` logic handles the
  remaining 5-node subset correctly without needing a config change.
- **`config/L2_olink_v2.R`** — not edited, not even reviewed for
  Olink-specific redundancy (e.g. related biomarker isoforms). L2 was never
  affected by anything in this session beyond the two shared `.Rmd` fixes
  (threshold logic, WBC blacklist — which correctly no-ops for L2 since it
  has none of the relevant CBC nodes).
- **`bn_learning/` (bootstrap pipeline)** — entirely untouched. The
  deterministic-column exclusions made to L1a's config in this session are
  **not** mirrored in the bootstrap pipeline's version of the same config,
  if one exists separately. If bootstrap-vs-CV comparisons are planned, this
  is a real inconsistency to resolve first.
- **`chloride`/`sodium`/`potassium` cycle warnings** (electrolyte cluster,
  correlation-driven but not deterministic) — deliberately left unfixed.
  Judged to be real (if noisy) physiological co-variation rather than a
  calculated-metric artifact; blacklisting would suppress potentially
  genuine signal. `averaged.network()`'s built-in cycle-skip handles these
  gracefully as-is.
- **`neutrophils -> sodium` cycle warning** (L1a) — checked against raw
  data (correlation only -0.365, no data-quality red flags), concluded to
  be sampling noise at n=81 rather than a systematic issue. Left as-is for
  the same reason as the electrolyte cluster.
- **L2's CV loss variance** (sd/mean ≈ 20%, higher than L1a/L1b's ~5-7%
  post-fix) — flagged as worth a closer look (checked, appears to be a
  smooth spread across folds rather than a single dominant outlier fold,
  unlike L1a's original problem) but no action taken.

---

## Open items / stale inconsistencies noticed but not resolved

- **L1a config header comment mismatch:** the file's top-of-file comment
  says "3 roots + 47 clinical continuous nodes," and an inline comment
  (pre-existing, before this session's edits) said "47 = 53 - 6
  exclude_candidates" — but `CONTINUOUS_MAP` only ever had 44 entries
  before this session's changes (now 34). This discrepancy predates this
  session and was not investigated further.
