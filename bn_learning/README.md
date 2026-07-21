# bn_learning Pipeline — Bootstrap Structure Learning Results

## Overview

CLG (Conditional Linear Gaussian) Bayesian networks learned from the Sound Life
Cohort baseline multi-omics data (n=92). Three discrete root nodes
(age group, biological sex, CMV serostatus) drive continuous clinical and
molecular features. Structure learning uses bootstrap model averaging
(`boot.strength()`, 1000 replicates with tabu search) rather than
cross-validation (see `bn_cv/` for the CV pipeline).

## What each network is

| Network | Modality | Description | Config |
|---|---|---|---|
| **L1a** | Clinical | Full clinical panel — every clinical lab/measurement column not excluded for missingness or determinism | `config/L1a_clinical_full.R` |
| **L1a_v3** | Clinical | L1a with BMI/height excluded (missingness) + 10 deterministic/derived columns removed (perc_*, mch/mchc/mcv, non_hdl, hdl_ratio) | `config/L1a_clinical_full_v3.R` |
| **L1b** | Clinical | Curated clinical panel — literature-selected ~22 clinically meaningful markers | `config/L1b_clinical_curated.R` |
| **L2** | Olink | Proteomics — 35-protein inflammation/immune panel via Olink NPX assay | `config/L2_olink.R` |
| **L2_v2** | Olink | Olink with MMP9 excluded (high missingness; n rises 72→88) | `config/L2_olink_v2.R` |
| **L3a** | Combined | L1a (full clinical) + L2 (Olink) merged into one joint network | `config/L3a_clinical_full_olink.R` |
| **L3b** | Combined | L1b (curated clinical) + L2 (Olink) merged into one joint network | `config/L3b_clinical_curated_olink.R` |
| **L4** | Whole blood | Bulk RNA-seq scored against 50 MSigDB Hallmark gene sets — pathway-level activity summaries | `config/L4_wb_pathways.R` |
| **L5** | Immunophenotyping | Flow/mass cytometry cell-type frequencies, 3-level gating hierarchy with ALR transform | `config/L5_freq_immunophenotype.R` |
| **L6** | PBMC signaling | Phospho-flow cytometry — PBMCs stimulated with 6 signals across 24 cell types (144 nodes) | `config/L6_pb_signaling.R` |
| **L_all** | Kitchen sink | All modalities combined, programmatic column discovery, maxp=5 | `config/L_all.R` |

All networks share the same 3 root nodes: `age_group`, `sex`, `cmv` (CMV serostatus).

---

## Cross-cutting fixes (ported from bn_cv, 2026-07-20)

| Fix | What it does | Applies to |
|---|---|---|
| Empirical-valley conservative threshold | Replaced fixed 0.85 cutoff with per-network histogram-based valley in bootstrap edge-strength distribution (0.3–0.95 range, bin width 0.05) | All networks (re-rendered after fix) |
| WBC-cluster blacklist (Rule 3) | Forbids direct edges between WBC and its differential white-cell-count components via `intersect()` | L1a, L1b (partial), L3a, L3b; no-ops elsewhere |
| Hierarchy blacklist (Rule 4) | Forbids direct parent-child edges using config-supplied `PARENT_MAP` | L5, L_all (freq columns only); no-ops for all others |
| `maxp` parameter | Caps parent count per node in `hc()`/`tabu()`/`mmhc()`/`boot.strength()` — prevents near-singular local fits | All networks (default 10; L_all uses 5) |
| `cache=TRUE` removed | knitr cache keys on chunk code not `params$config` — silently reused first config's bootstrap results for all subsequent networks | All notebooks |

---

## L1a — Clinical (full panel, original)

**49 nodes** (3 roots + 46 continuous), **69/92 samples**.

**Status:** Original config with all derived columns still included. Superseded
by L1a_v3 — kept for historical comparison only.

**Results:**
```
Tabu: 138 edges
avg_50 (t=0.50): 88 edges
Edges with strength >= 0.85: 33 (undirected)
```

**Dominated by trivial/derived edges:** bmi↔weight (1.000),
basophils↔perc_basophils (1.000), eosinophils↔perc_eosinophils (1.000),
mch↔mchc (1.000), mch↔mcv (1.000), ldl↔non_hdl (1.000),
perc_lymphocytes↔perc_neutrophils (1.000). These are exact arithmetic
functions of other columns, not biology.

---

## L1a_v3 — Clinical (full panel, cleaned)

**37 nodes** (3 roots + 34 continuous), **81/92 samples**.

**Fixes applied:** BMI/height excluded (missingness), 10 deterministic columns
removed (mch/mchc/mcv, 5 perc_*, non_hdl, hdl_ratio), WBC-cluster blacklist,
empirical valley threshold, maxp=10.

**Results:**
```
hc: 80 edges   tabu: 80 edges   mmhc: 29 edges
avg_opt: 55 edges
avg_50 (t=0.50): 55 edges
maxp: 10
Edges with strength >= 0.85: 18 (undirected)
```

**Top edges (strength >= 0.85):**
| Edge | Strength | Biology |
|---|---|---|
| hematocrit — hemoglobin | 1.000 | Real co-regulation (both RBC-dependent) |
| alt — ast | 1.000 | Liver transaminases, biologically coupled |
| ldl — total_chol | 1.000 | Friedewald-related but not exactly deterministic |
| mpv — platelets | 0.996 | Inverse relationship (thrombopoiesis) |
| hematocrit — rbc | 0.994 | RBC mass determines hematocrit |
| hdl — triglycerides | 0.990 | Metabolic syndrome axis |
| globulin — protein | 0.986 | protein ≈ albumin + globulin (near-deterministic) |
| albumin — calcium | 0.983 | ~45% serum calcium is albumin-bound |
| rdw — rbc | 0.976 | RBC production/turnover |
| sex → crp | 0.876 | Sex differences in inflammatory markers |

Much cleaner than L1a — derived-column edges eliminated. Still contains some
near-deterministic pairs (globulin↔protein, ldl↔total_chol) that could be
excluded in a future iteration.

---

## L1b — Clinical (curated panel)

**25 nodes** (3 roots + 22 continuous), **87/92 samples**.

**Fixes applied:** WBC-cluster blacklist (partial — no basophils in this panel),
empirical valley threshold, maxp=10.

**Results:**
```
hc: 49 edges   tabu: 49 edges   mmhc: 15 edges
avg_opt: 36 edges
avg_50 (t=0.50): 36 edges
maxp: 10
Edges with strength >= 0.85: 7 (undirected)
```

**Top edges (all with strength >= 0.85):**
| Edge | Strength | Biology |
|---|---|---|
| hemoglobin — hematocrit | 1.000 | RBC biology |
| alt — ast | 1.000 | Liver transaminases |
| hdl — triglycerides | 0.999 | Metabolic syndrome axis |
| albumin — calcium | 0.996 | Albumin-bound calcium |
| creatinine — bun | 0.931 | Renal function axis |
| esr — platelets | 0.860 | Both acute-phase markers |
| albumin — bun | 0.851 | Hepatorenal axis |

The cleanest, most interpretable network. Small node count (25) and high
sample efficiency (87/92) yield a sparse, stable graph. All high-strength
edges are well-established clinical relationships. CMV has no strong edges
here — consistent with CMV driving proteins/cells (L2), not clinical labs.

---

## L2 — Olink proteomics (original)

**38 nodes** (3 roots + 35 continuous), **72/92 samples**.

**Results (pre-stability-fix, old format):**
```
hc: 86 edges   tabu: 90 edges   mmhc: 37 edges
avg_50 (t=0.50): 61 edges
Edges with strength >= 0.85: 13 (undirected)
```

Superseded by L2_v2 with MMP9 excluded and stability fixes applied.

---

## L2_v2 — Olink proteomics (cleaned)

**37 nodes** (3 roots + 34 continuous), **88/92 samples**.

**Fixes applied:** MMP9 excluded (missingness), empirical valley threshold,
maxp=10.

**Results:**
```
hc: 89 edges   tabu: 91 edges   mmhc: 36 edges
avg_opt: 58 edges
avg_50 (t=0.50): 58 edges
maxp: 10
Edges with strength >= 0.85: 14 (undirected)
```

**Top edges (strength >= 0.85):**
| Edge | Strength | Biology |
|---|---|---|
| IL7 — VEGFA | 0.995 | Growth factor co-regulation |
| age → GDF15 | 0.991 | THE aging biomarker (mitokine) |
| sex → LEP | 0.991 | Leptin sex dimorphism |
| CXCL9 — CXCL10 | 0.990 | IFN-γ-induced chemokine pair |
| TNFRSF1B — HAVCR2 | 0.983 | TNF receptor / Tim-3 checkpoint |
| GZMA — LAG3 | 0.958 | Cytotoxic/checkpoint cluster |
| IL1RN — LEP | 0.947 | Anti-inflammatory / adipokine link |
| TNFRSF1A — TNFRSF1B | 0.936 | TNF receptor pair |
| CD8A — GZMA | 0.925 | CD8 T-cell cytotoxicity |
| CXCL10 — CXCL11 | 0.924 | IFN-γ chemokine cascade |
| TNFRSF1A — HAVCR2 | 0.903 | TNF/checkpoint subnetwork |
| IL6 — IL1RN | 0.888 | Pro/anti-inflammatory balance |
| CXCL9 — PDCD1 | 0.882 | Chemokine/checkpoint (PD-1) |
| IL18 — IL1RN | 0.872 | IL-1 family cross-regulation |

Rich immune biology recovered from data alone: IFN-γ chemokine cascade
(CXCL9→CXCL10→CXCL11), cytotoxic/checkpoint cluster (CD8A↔GZMA→LAG3),
TNF receptor subnetwork (TNFRSF1A↔TNFRSF1B↔HAVCR2). v2 confirms v1
structure: core edges stable despite +16 samples.

---

## L3a — Combined: Clinical (full) + Olink

**71 nodes** (3 roots + 34 clinical + 34 Olink), **81/92 samples** (ratio 0.88).

**Results:**
```
hc: 224 edges   tabu: 231 edges   mmhc: 70 edges
avg_opt: 104 edges
avg_50 (t=0.50): 104 edges
maxp: 10
Edges with strength >= 0.85: 23 (undirected)
```

**Issue:** Still contains the full clinical panel with near-deterministic pairs
(ldl↔total_chol at 1.000, globulin↔protein at 0.960, hematocrit↔rbc at 0.994)
that crowd the top of the strength ranking. The same derived-column fix applied
to L1a_v3 should be applied here for a clean comparison.

**Cross-modality edges (strength >= 0.50):**
| Edge | Strength | Biology |
|---|---|---|
| weight → LEP | 0.977 | Leptin from adipose tissue — direct biological link |
| protein ↔ MMP12 | 0.730 | Matrix metalloproteinase / serum protein |
| bun ↔ FGF21 | 0.715 | FGF21 is a liver/kidney metabolic hormone |
| alk_phos ↔ IL1RN | 0.672 | Bone/liver enzyme + anti-inflammatory cytokine |
| eosinophils ↔ LEP | 0.625 | Adipokine-eosinophil axis |
| CRP ↔ IL6 | 0.570 | Classic inflammation circuit (IL-6 drives hepatic CRP) |
| triglycerides ↔ VEGFA | 0.554 | Metabolic-vascular link |
| lymphocytes ↔ IL15 | 0.549 | IL-15 is key lymphocyte survival cytokine |
| glucose ↔ IL10 | 0.543 | Anti-inflammatory / metabolic link |

The CRP↔IL6 edge (0.570, direction ~50/50) is notable — the question "does
IL6 mediate CRP?" is now testable in this combined network, though direction
is not resolvable with this sample size.

---

## L3b — Combined: Clinical (curated) + Olink

**59 nodes** (3 roots + 22 clinical + 34 Olink), **87/92 samples** (ratio 0.68).

**Results:**
```
hc: 179 edges   tabu: 182 edges   mmhc: 54 edges
avg_opt: 89 edges
avg_50 (t=0.50): 89 edges
maxp: 10
Edges with strength >= 0.85: 15 (undirected)
```

**The best combined network.** No trivial derived edges in the top hits —
the curated clinical panel already excluded those. More samples (87 vs 81)
than L3a. All high-strength edges are biologically real:

**Top edges (strength >= 0.85):**
| Edge | Strength | Biology |
|---|---|---|
| hemoglobin — hematocrit | 1.000 | RBC biology |
| alt — ast | 1.000 | Liver transaminases |
| albumin — calcium | 0.995 | Albumin-bound calcium |
| IL7 — VEGFA | 0.991 | Growth factor co-regulation |
| hdl — triglycerides | 0.990 | Metabolic syndrome axis |
| CXCL9 — CXCL10 | 0.986 | IFN-γ chemokines |
| age → GDF15 | 0.979 | Aging biomarker |
| TNFRSF1B — HAVCR2 | 0.965 | TNF/checkpoint |
| sex → LEP | 0.964 | Leptin sex dimorphism |
| GZMA — LAG3 | 0.939 | Cytotoxic/checkpoint |
| CD8A — GZMA | 0.937 | CD8 T-cell cytotoxicity |
| IL1RN — LEP | 0.904 | Anti-inflammatory/adipokine |
| TNFRSF1A — TNFRSF1B | 0.890 | TNF receptor pair |
| TNFRSF1A — HAVCR2 | 0.887 | TNF/checkpoint subnetwork |
| CXCL10 — CXCL11 | 0.882 | IFN-γ cascade |

**Cross-modality edges (strength >= 0.50):**
| Edge | Strength | Biology |
|---|---|---|
| bun ↔ FGF21 | 0.700 | Liver/kidney metabolic hormone |
| CRP ↔ VEGFA | 0.675 | Inflammation-vascular link |
| eosinophils ↔ LEP | 0.672 | Adipokine-eosinophil axis |
| triglycerides ↔ VEGFA | 0.656 | Metabolic-vascular |
| neutrophils ↔ CXCL8 | 0.637 | IL-8 is primary neutrophil chemoattractant |
| glucose ↔ IL10 | 0.624 | Anti-inflammatory/metabolic |
| triglycerides ↔ IGFBP7 | 0.618 | Metabolic signaling |
| CRP ↔ IL6 | 0.599 | Classic inflammation circuit |
| creatinine ↔ IL18 | 0.590 | IL-18 cleared renally |
| alt ↔ IL18 | 0.580 | Liver damage / inflammasome |

The combined network successfully recovers both within-modality structure
(clinical pairs from L1b, Olink clusters from L2) and cross-modality
connections. neutrophils↔CXCL8 and creatinine↔IL18 are particularly clean
biological findings.

---

## L4 — Whole blood (Hallmark pathway scores)

**53 nodes** (3 roots + 50 continuous), **86/92 samples** (ratio 0.62).

**Fixes applied:** `wb.low_lib_flag` excluded (boolean QC flag, not a pathway
score). No deterministic-collinearity issues — Hallmark pathway scores can
share member genes and correlate biologically, but aren't exact arithmetic
functions of one another.

**Results:**
```
hc: 139 edges   tabu: 139 edges   mmhc: 46 edges
avg_50 (t=0.50): 86 edges
maxp: 10
Edges with strength >= 0.85: 13 (undirected)
```

**Top edges (strength >= 0.85):**
| Edge | Strength | Biology |
|---|---|---|
| e2f_targets — g2_m_checkpoint | 1.000 | Both cell-cycle regulation pathways |
| interferon_alpha — interferon_gamma | 1.000 | Shared IFN signaling / STAT1 convergence |
| estrogen_response_early — estrogen_response_late | 1.000 | Same hormone, temporal phases |
| g2_m_checkpoint — mitotic_spindle | 0.997 | Cell division machinery |
| uv_response_dn — epithelial_mesenchymal_transition | 0.964 | Stress response / tissue remodeling |
| epithelial_mesenchymal_transition — apical_junction | 0.964 | Cell structure / tissue architecture |
| complement — coagulation | 0.963 | Both serine protease cascades, shared components |
| hypoxia — protein_secretion | 0.958 | Hypoxia drives UPR/secretory pathway |
| inflammatory_response — il_6_jak_stat3 | 0.901 | IL-6/JAK/STAT3 is a core inflammatory signaling axis |
| cholesterol_homeostasis — apoptosis | 0.892 | Lipid metabolism / cell death crosstalk |
| e2f_targets — myc_targets_v1 | 0.890 | MYC drives E2F transcription |
| inflammatory_response — kras_signaling_up | 0.869 | KRAS activates NF-κB / inflammatory programs |
| mitotic_spindle — adipogenesis | 0.862 | Cell division in differentiating adipocytes |

**Root edges:** Nearly no root-driven signal — only sex → reactive_oxygen_species
(0.519). Age and CMV have no strong pathway-level effects in whole blood,
suggesting their influence operates at the cell-type or protein level
(L2/L5), not via bulk transcriptomic pathway scores.

The pathway-level network is dominated by biologically coherent functional
modules: a cell-cycle cluster (e2f/g2m/mitotic_spindle/myc), an immune
signaling cluster (IFN-α/γ, inflammatory_response, IL-6/JAK/STAT3,
complement/coagulation), and a tissue-remodeling cluster (EMT, apical
junction, UV response).

---

## L5 — Immune cell frequencies (3-level gating hierarchy)

**91 nodes** (3 roots + 88 continuous after ALR transform), **88/92 samples**
(ratio 1.03). Original config has 112 nodes (109 continuous); the ALR
transform reduces sibling groups and drops reference nodes.

**Fixes applied:**
- 12 exact-duplicate columns removed (`.2`/`.3` suffixed variants, source-data
  concatenation artifact — same fix as bn_cv).
- ALR transform for compositional hierarchy: for each sibling group under a
  shared parent, all but the highest-mean-abundance sibling replaced with
  `log(sibling) - log(reference)`, reference dropped. Extended to treat the
  9 L1 categories as a top-level sibling group.
- Hierarchy blacklist (Rule 4) via `PARENT_MAP` — forbids direct parent-child
  edges in the gating tree.
- `maxp=10`.

**Results:**
```
hc: 428 edges   tabu: 431 edges   mmhc: 90 edges
avg_50 (t=0.50): 153 edges
maxp: 10
Edges with strength >= 0.85: 52 (undirected)
```

**Structural note:** Many of the strongest edges (str=1.000) connect nodes
across adjacent gating-tree levels where only one child exists under a parent
(e.g. `l2_erythrocyte — l1_erythrocyte_alr`, `l2_platelet — l1_platelet_alr`).
These are not trivial/derived — the ALR transform means they measure different
things (absolute vs. relative abundance) — but they reflect the gating
hierarchy's structure more than independent biology. The interesting edges are
the within-level connections:

**Notable within-level edges (strength >= 0.85):**
| Edge | Strength | Biology |
|---|---|---|
| l2_effector_b_cell_alr — l2_memory_b_cell_alr | 1.000 | B-cell differentiation axis |
| l3_isg_cd14_monocyte_alr — l3_isg_cd16_monocyte_alr | 1.000 | Interferon-stimulated monocyte subtypes |
| l3_naive_vd1_gdt_alr — l3_sox4_vd1_gdt_alr | 1.000 | Vδ1 γδ T-cell maturation states |
| l1_dc_alr — l1_monocyte_alr | 0.999 | Shared myeloid progenitor |
| l3_sox4_naive_cd4 — l3_sox4_naive_cd8 | 0.996 | SOX4+ naive T-cell co-regulation |
| l2_memory_cd4_t_cell_alr — l2_treg_alr | 0.995 | Tregs arise from memory CD4 pool |

**Root-connected edges (strength >= 0.5):**
| Edge | Strength | Biology |
|---|---|---|
| cmv → klrf1_gzmb_cd27_memory_cd4_t_cell_alr | 1.000 | CMV drives cytotoxic memory CD4 expansion — strongest root effect in any network |
| age → sox4_naive_cd4_t_cell_alr | 0.816 | Thymic involution reduces naive T-cells with age |
| cmv → adaptive_nk_cell_alr | 0.800 | CMV induces adaptive/memory NK cells |
| sex → early_memory_b_cell_alr | 0.702 | Sex differences in B-cell maturation |
| age → cd16_monocyte_alr | 0.661 | Age-related monocyte subset shift |
| cmv → klrf1_effector_vd1_gdt_alt_alr | 0.550 | CMV expands Vδ1 γδ T effectors |
| age → gdt_alr | 0.508 | Age-related γδ T-cell decline |

**CMV dominates the root effects** in L5 — 3 of the 7 root edges are
CMV-driven, consistent with CMV's known large-scale reshaping of the
cellular immune compartment (particularly memory CD4, adaptive NK, and
γδ T cells). This is the modality where CMV has its strongest signal,
complementing the protein-level CMV effects seen in L2 (KLRD1, TNF).

---

## L6 — PBMC functional/signaling response

**147 nodes** (3 roots + 144 continuous = 24 cell types × 6 stimulation
conditions), **77/92 samples** (ratio 1.91 — highest of any network).

**Fixes applied:** Cell-type selection restricted to 24 of 65 measured types
(5% per-cell-type missingness cutoff keeps n at 77/92). 65 `*_low_ncells_flag`
QC columns excluded. maxp=10.

**Results:**
```
hc: 719 edges   tabu: 721 edges   mmhc: 141 edges
avg_50 (t=0.50): 216 edges
maxp: 10
Edges with strength >= 0.85: 60 (undirected)
```

**Dominant pattern: IFN-α ↔ IFN-γ within-cell-type edges.** ~20 of the 24
cell types show their IFN-α and IFN-γ conditions linked at strength 1.000
or near it. This reflects both interferon types converging on shared STAT1
phosphorylation — real biology, not an artifact, but it means the strongest
edges are structurally predictable rather than novel. The same pattern was
confirmed in bn_cv and intentionally kept (max correlation 0.976, below the
exact-duplicate threshold).

**Same cell type, different conditions (top edges):**
| Edge | Strength | Pattern |
|---|---|---|
| ~20 cell types: *_ifn_alpha — *_ifn_gamma | 1.000 | Shared STAT1 signaling |
| pdc_il2_stat5 — pdc_il6_jak_stat3 | 0.994 | Cross-STAT signaling in pDCs |

**Different cell types, same condition (more biologically interesting):**
| Edge | Strength | Biology |
|---|---|---|
| core_cd14_monocyte — isg_cd14_monocyte (TNFα/NF-κB) | 1.000 | Monocyte subtype co-activation |
| core_cd14_monocyte — isg_cd14_monocyte (IL-6/JAK/STAT3) | 0.999 | Same monocyte pair, different pathway |
| core_cd14_monocyte — isg_cd14_monocyte (IL-2/STAT5) | 0.999 | Same monocyte pair, third pathway |
| core_naive_cd8 — core_naive_cd4 (IFN-γ) | 0.996 | Naive T-cell co-regulation |
| core_cd14_monocyte — isg_cd14_monocyte (inflammatory) | 0.996 | Same monocyte pair, fourth pathway |
| core_naive_cd8 — core_naive_cd4 (TNFα/NF-κB) | 0.994 | Naive T-cell co-activation |
| core_naive_cd8 — core_naive_cd4 (IFN-α) | 0.983 | Naive T-cells, third condition |
| cm_cd4 — gzmb_cd27_em_cd4 (TNFα/NF-κB) | 0.981 | CD4 memory/effector axis |
| core_naive_cd8 — core_naive_cd4 (IL-2/STAT5) | 0.971 | Naive T-cells, fourth condition |
| cd14_cdc2 — hla_drhi_cdc2 (TNFα/NF-κB) | 0.952 | DC subtype co-activation |
| core_naive_b — core_memory_b (TNFα/NF-κB) | 0.922 | B-cell maturation axis |
| gzmk_cd56dim_nk — cd8_mait (TNFα/NF-κB) | 0.915 | Innate-like lymphocyte co-response |

**Key structural insight:** The core_cd14_monocyte ↔ isg_cd14_monocyte pair
is linked across ALL tested conditions (str ≥ 0.996 for TNFα/NF-κB,
IL-6/JAK/STAT3, IL-2/STAT5, inflammatory). These two monocyte subtypes
respond in lockstep regardless of stimulus, suggesting they share the same
upstream activation pathway or are subsets of the same functional state.

**Root edges (strength >= 0.5):** Only 2, both age-related:
| Edge | Strength | Biology |
|---|---|---|
| age → core_naive_cd8_il2_stat5 | 0.654 | Age reduces naive CD8 IL-2 responsiveness |
| age → core_naive_cd8_il6_jak_stat3 | 0.548 | Age reduces naive CD8 STAT3 signaling |

CMV and sex have no strong effects at the signaling level — CMV reshapes the
cellular compartment (L5) and protein milieu (L2) but does not detectably
alter per-cell signaling responses. Age's effect is limited to naive CD8
T-cells, consistent with thymic involution reducing naive T-cell functional
capacity.

**Runtime note:** Bootstrap took ~3 hours (1000 replicates × tabu on 147
nodes). Future reruns could use `maxp=5` and/or `R=200` to reduce runtime
without substantially affecting high-strength edge detection.

---

## L_all — Kitchen sink (all modalities)

**352 nodes** (3 roots + 349 continuous), **75/92 samples** (ratio 4.69 —
by far the highest of any network). Programmatic column discovery from
`bn_ready_baseline.csv` using `FEATURE_PREFIXES`, 5% per-column missingness
threshold before complete-case, `maxp=5`.

**Modality breakdown:** 34 Olink + 88 freq (ALR-transformed) + 145 pb
signaling + clinical (aliased with `bc_`, `chem_`, `lip_`, `infl_` prefixes).
50 whole blood pathway scores also included.

**Results:**
```
hc: 1533 edges   tabu: 1537 edges   mmhc: 350 edges
avg_opt: 305 edges
avg_50 (t=0.50): 305 edges
maxp: 5
Edges with strength >= 0.85: 118 (undirected)
```

**High-strength edges recapitulate individual-network findings.** The
strongest edges are the same ones discovered in the single-modality networks:
- Clinical: hematocrit↔hemoglobin, neutrophils↔wbc, alt↔ast, ldl↔total_chol (all 1.000)
- Freq: gating-tree hierarchy edges (l2↔l1_alr, l3↔l2_alr at 1.000)
- PB: IFN-α↔IFN-γ within-cell-type pairs, core↔isg monocyte co-activation (all 1.000)

**Root-connected edges (strength >= 0.5):**
| Edge | Strength | Modality | Biology |
|---|---|---|---|
| cmv → l3_klrf1_gzmb_cd27_memory_cd4_t_cell_alr | 0.998 | freq | CMV cytotoxic memory CD4 expansion |
| sex → LEP | 0.893 | olink | Leptin sex dimorphism |
| sex → chem_creatinine | 0.881 | clinical | Sex differences in muscle mass/creatinine |
| age → GDF15 | 0.824 | olink | THE aging biomarker |
| age → l3_sox4_naive_cd4_t_cell_alr | 0.604 | freq | Thymic involution |
| cmv → l3_gzmk_cd27_em_cd8_t_cell_alt_alr | 0.586 | freq | CMV CD8 effector expansion |
| cmv → l3_adaptive_nk_cell_alr | 0.514 | freq | CMV adaptive NK cells |
| sex → IGFBP7 | 0.508 | olink | Sex-linked insulin/growth factor binding |

**Key finding:** The root effects survive the kitchen-sink model and remain
modality-specific — CMV drives freq (cellular composition), sex/age drive
olink (protein levels) and clinical labs. No novel cross-modality root effects
emerge at this scale, validating the single-modality network results.

**Cross-modality edges are sparse** at the conservative threshold. The
within-modality structure from L1b–L6 dominates. This is expected at a 4.7:1
node-to-sample ratio with maxp=5 — the aggressive parent cap correctly
prevents overfitting but also limits the network's ability to discover weaker
cross-modality connections. The individual L3a/L3b combined networks remain
better suited for cross-modality investigation.

**Runtime note:** Bootstrap took several hours (1000 replicates × tabu on 352
nodes with maxp=5).

---

## Summary table

| Network | Nodes | Samples | Ratio | Tabu edges | avg_opt edges | Str≥0.85 | Status |
|---|---|---|---|---|---|---|---|
| L1a | 49 | 69 | 0.71 | 138 | — | 33 | Superseded by v3 |
| L1a_v3 | 37 | 81 | 0.46 | 80 | 55 | 18 | ✅ Clean |
| L1b | 25 | 87 | 0.29 | 49 | 36 | 7 | ✅ Cleanest |
| L2 | 38 | 72 | 0.53 | 90 | — | 13 | Superseded by v2 |
| L2_v2 | 37 | 88 | 0.42 | 91 | 58 | 14 | ✅ Clean |
| L3a | 71 | 81 | 0.88 | 231 | 104 | 23 | ⚠️ Needs derived-column fix |
| L3b | 59 | 87 | 0.68 | 182 | 89 | 15 | ✅ Best combined |
| L4 | 53 | 86 | 0.62 | 139 | 86 | 13 | ✅ Clean |
| L5 | 91 | 88 | 1.03 | 431 | 153 | 52 | ✅ Clean |
| L6 | 147 | 77 | 1.91 | 721 | 216 | 60 | ✅ Clean |
| L_all | 352 | 75 | 4.69 | 1537 | 305 | 118 | ✅ Clean (maxp=5) |

---

## Directory layout

```
bn_learning/
├── 00_install_packages.R           # One-time R package setup
├── 01_data_prep.Rmd                # Standard data prep (L1a/L1b/L2/L3a/L3b/L4/L6)
├── 01_data_prep_L5.Rmd             # L5-specific: ALR transform for freq hierarchy
├── 01_data_prep_all.Rmd            # L_all-specific: programmatic column discovery
├── 02_structure_learning.Rmd       # Shared: hc/tabu/mmhc + bootstrap (parameterized)
├── 03_model_comparison.Rmd         # Cross-level comparison
├── export_edge_lists.R             # Standalone: extract edge CSVs from _models.rds
├── export_network_plots.R          # Standalone: export DAG PNGs from _models.rds
├── config/                         # Network-level configs (sourced by notebooks)
├── output/
│   ├── figures/                    # Edge strength histograms
│   ├── mappings/                   # *_node_mapping.csv (alias <-> csv column)
│   ├── models/                     # .rds model + data + boot_str files
│   └── networks/
│       ├── csv_networks/           # *_literature_edges.csv (full edge lists)
│       └── visualized_networks/    # *_tabu.png, *_avg_opt.png, *_avg_conservative.png
└── rendered/                       # HTML outputs per network level
    ├── L1a/                        # v1, v2, v3 renders
    ├── L1b/
    ├── L2/
    ├── L3a/
    ├── L3b/
    ├── L4/
    ├── L5/
    ├── L6/
    ├── L_all/
    └── comparison/
```

## Notebooks

| Notebook | Params | Purpose |
|----------|--------|---------|
| `01_data_prep.Rmd` | `config` | Load data, backfill CMV, transform, subset, complete-case |
| `01_data_prep_L5.Rmd` | `config` | Same + ALR transform for compositional freq hierarchy |
| `01_data_prep_all.Rmd` | `config` | Programmatic column discovery + missingness filter |
| `02_structure_learning.Rmd` | `config` | hc/tabu/mmhc + 1000-replicate bootstrap model averaging |
| `03_model_comparison.Rmd` | `version` | Compare algorithms, cross-level analysis |

## Algorithms

- `bn_hc`: hill-climbing, `score="bic-cg"` (greedy)
- `bn_tabu`: tabu search, `score="bic-cg"` (primary — escapes local optima)
- `bn_mmhc`: hybrid (constraint skeleton + score-based orientation)
- `boot_str`: bootstrap edge strengths (1000 replicates, tabu-based)
- `avg_opt`: averaged network at optimal threshold (bnlearn heuristic)
- `avg_conservative`: averaged network at empirical valley threshold
- `avg_50`: averaged network at t=0.50

## Key design decisions

- **CLG**: discrete roots cannot have continuous parents (enforced by blacklist)
- **CMV backfilled**: only 9/92 measured at baseline; 88/92 have CMV from other visits (stable trait); 4 subjects (BR1049, BR1056, BR1058, BR1059) have no CMV at any visit and are dropped at complete-case
- **Log-transformed**: CRP, ESR, ALT, AST, triglycerides, glucose (clinical networks only)
- **Olink NPX**: already log2-scale, no additional transform
- **Hallmark pathway scores**: already normalized enrichment scores, no additional transform
- **Complete-case**: rows with any NA dropped (varies by network level)
- **NLR excluded**: deterministic function of neutrophil + lymphocyte counts (all clinical networks)
- **10 derived columns excluded** (L1a_v3, L3a DETERMINISTIC_EXCLUDE): mch/mchc/mcv, 5 perc_*, non_hdl, hdl_ratio

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
