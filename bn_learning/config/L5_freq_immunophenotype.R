# config/L5_freq_immunophenotype.R
# L5: Immune cell frequency-only — 3 roots + 109 cell-type frequency nodes
# Flow/mass cytometry immunophenotyping: for each subject, the percentage of
# cells falling into a defined cell-type bucket, at three nested levels of
# granularity (L1: 9 broad lineages, L2: 29 finer subtypes, L3: originally
# 83 finest subtypes).
#
# STRUCTURAL NOTE: this modality is fundamentally different from clinical/
# Olink/wb. L1/L2/L3 form a gating hierarchy, not independent measurements
# -- L2 subtypes are subsets of their L1 parent lineage, L3 subsets of their
# L2 parent. This creates the same kind of collinearity risk as WBC vs. its
# differential counts in the clinical modality, but systemically across the
# whole node set rather than one isolated cluster. Handled here via
# PARENT_MAP (below) + a corresponding blacklist rule in
# 02_structure_learning.Rmd (forbids direct parent-child edges; does NOT
# forbid edges between siblings or non-adjacent levels, since those aren't
# guaranteed redundant the way a direct parent-child pair is).

NETWORK_ID <- "L5_freq_immunophenotype"
SEED <- 46L

EXCLUDE_CANDIDATES <- character(0)  # 0% missingness across all 121 raw columns

# DETERMINISTIC_EXCLUDE: 12 columns removed as EXACT duplicates, not derived
# formulas (a stronger case than NLR/MCHC -- these are byte-identical copies,
# confirmed against bn_ready_baseline.csv). Six L3 base names each had 4
# suffixed variants (base, .1, .2, .3) from what appears to be a source-data
# concatenation issue: base==.2 exactly and .1==.3 exactly, in all 92
# subjects, across all 6 groups, with no exceptions. Kept `base` and `.1`
# (renamed to *_alt below, since what actually distinguishes them from
# `base` is unknown -- verify against the original data source if this
# matters for interpretation) and dropped `.2`/`.3` as pure duplicates.
DETERMINISTIC_EXCLUDE <- c(
  "freq.L3_cd27_effector_b_cell.2",              "freq.L3_cd27_effector_b_cell.3",
  "freq.L3_gzmb_cd27_em_cd4_t_cell.2",            "freq.L3_gzmb_cd27_em_cd4_t_cell.3",
  "freq.L3_gzmk_cd27_em_cd8_t_cell.2",            "freq.L3_gzmk_cd27_em_cd8_t_cell.3",
  "freq.L3_gzmk_cd56dim_nk_cell.2",               "freq.L3_gzmk_cd56dim_nk_cell.3",
  "freq.L3_klrf1_effector_vd1_gdt.2",             "freq.L3_klrf1_effector_vd1_gdt.3",
  "freq.L3_klrf1_gzmb_cd27_em_cd8_t_cell.2",      "freq.L3_klrf1_gzmb_cd27_em_cd8_t_cell.3"
)

# KNOWN GAP: l3_baeomap_cell (below) is real, non-duplicate data, but is not
# a standard immunology term I can confidently map to an L2 parent -- it is
# included as a node but deliberately OMITTED from PARENT_MAP, meaning it
# gets no blacklist protection. If it causes cycle warnings later, that's
# expected and lower-priority than fabricating a wrong biological claim.
# Verify this term against whoever generated the freq modality before
# treating its position in the hierarchy as resolved.

ROOT_MAP <- c(
  age_group = "subject.ageGroup",
  sex       = "subject.biologicalSex",
  cmv       = "cmv.igg_serology_interpretation"
)

# Continuous cell-frequency nodes (109 = 9 L1 + 29 L2 + 71 L3)
CONTINUOUS_MAP <- c(
  # --- L1: broad lineages (9) ---
  l1_b_cell               = "freq.L1_b_cell",
  l1_dc                   = "freq.L1_dc",
  l1_erythrocyte          = "freq.L1_erythrocyte",
  l1_ilc                  = "freq.L1_ilc",
  l1_monocyte             = "freq.L1_monocyte",
  l1_nk_cell              = "freq.L1_nk_cell",
  l1_platelet             = "freq.L1_platelet",
  l1_progenitor_cell      = "freq.L1_progenitor_cell",
  l1_t_cell               = "freq.L1_t_cell",
  # --- L2: finer subtypes (29) ---
  l2_asdc                 = "freq.L2_asdc",
  l2_cd14_monocyte        = "freq.L2_cd14_monocyte",
  l2_cd16_monocyte        = "freq.L2_cd16_monocyte",
  l2_cd56bright_nk_cell   = "freq.L2_cd56bright_nk_cell",
  l2_cd56dim_nk_cell      = "freq.L2_cd56dim_nk_cell",
  l2_cd8aa                = "freq.L2_cd8aa",
  l2_cdc1                 = "freq.L2_cdc1",
  l2_cdc2                 = "freq.L2_cdc2",
  l2_dn_t_cell            = "freq.L2_dn_t_cell",
  l2_effector_b_cell      = "freq.L2_effector_b_cell",
  l2_erythrocyte          = "freq.L2_erythrocyte",
  l2_gdt                  = "freq.L2_gdt",
  l2_ilc                  = "freq.L2_ilc",
  l2_intermediate_monocyte = "freq.L2_intermediate_monocyte",
  l2_mait                 = "freq.L2_mait",
  l2_memory_b_cell        = "freq.L2_memory_b_cell",
  l2_memory_cd4_t_cell    = "freq.L2_memory_cd4_t_cell",
  l2_memory_cd8_t_cell    = "freq.L2_memory_cd8_t_cell",
  l2_naive_b_cell         = "freq.L2_naive_b_cell",
  l2_naive_cd4_t_cell     = "freq.L2_naive_cd4_t_cell",
  l2_naive_cd8_t_cell     = "freq.L2_naive_cd8_t_cell",
  l2_pdc                  = "freq.L2_pdc",
  l2_plasma_cell          = "freq.L2_plasma_cell",
  l2_platelet             = "freq.L2_platelet",
  l2_progenitor_cell      = "freq.L2_progenitor_cell",
  l2_proliferating_nk_cell = "freq.L2_proliferating_nk_cell",
  l2_proliferating_t_cell = "freq.L2_proliferating_t_cell",
  l2_transitional_b_cell  = "freq.L2_transitional_b_cell",
  l2_treg                 = "freq.L2_treg",
  # --- L3: finest subtypes (70, deduplicated from 83; baeomap_cell unmapped) ---
  l3_activated_memory_b_cell             = "freq.L3_activated_memory_b_cell",
  l3_asdc                                = "freq.L3_asdc",
  l3_c1q_cd16_monocyte                   = "freq.L3_c1q_cd16_monocyte",
  l3_cd14_cdc2                           = "freq.L3_cd14_cdc2",
  l3_cd27_effector_b_cell                = "freq.L3_cd27_effector_b_cell",
  l3_cd27_effector_b_cell_alt            = "freq.L3_cd27_effector_b_cell.1",
  l3_cd4_mait                            = "freq.L3_cd4_mait",
  l3_cd56bright_nk_cell                  = "freq.L3_cd56bright_nk_cell",
  l3_cd8_mait                            = "freq.L3_cd8_mait",
  l3_cd8aa                               = "freq.L3_cd8aa",
  l3_cd95_memory_b_cell                  = "freq.L3_cd95_memory_b_cell",
  l3_cdc1                                = "freq.L3_cdc1",
  l3_core_cd14_monocyte                  = "freq.L3_core_cd14_monocyte",
  l3_core_cd16_monocyte                  = "freq.L3_core_cd16_monocyte",
  l3_core_memory_b_cell                  = "freq.L3_core_memory_b_cell",
  l3_core_naive_b_cell                   = "freq.L3_core_naive_b_cell",
  l3_core_naive_cd4_t_cell               = "freq.L3_core_naive_cd4_t_cell",
  l3_core_naive_cd8_t_cell               = "freq.L3_core_naive_cd8_t_cell",
  l3_dn_t_cell                           = "freq.L3_dn_t_cell",
  l3_early_memory_b_cell                 = "freq.L3_early_memory_b_cell",
  l3_erythrocyte                         = "freq.L3_erythrocyte",
  l3_gzmb_vd2_gdt                        = "freq.L3_gzmb_vd2_gdt",
  l3_gzmk_cd56dim_nk_cell                = "freq.L3_gzmk_cd56dim_nk_cell",
  l3_gzmk_cd56dim_nk_cell_alt            = "freq.L3_gzmk_cd56dim_nk_cell.1",
  l3_gzmk_vd2_gdt                        = "freq.L3_gzmk_vd2_gdt",
  l3_hla_drhi_cdc2                       = "freq.L3_hla_drhi_cdc2",
  l3_il1b_cd14_monocyte                  = "freq.L3_il1b_cd14_monocyte",
  l3_ilc                                 = "freq.L3_ilc",
  l3_intermediate_monocyte               = "freq.L3_intermediate_monocyte",
  l3_isg_cd14_monocyte                   = "freq.L3_isg_cd14_monocyte",
  l3_isg_cd16_monocyte                   = "freq.L3_isg_cd16_monocyte",
  l3_isg_cd56dim_nk_cell                 = "freq.L3_isg_cd56dim_nk_cell",
  l3_isg_cdc2                            = "freq.L3_isg_cdc2",
  l3_isg_mait                            = "freq.L3_isg_mait",
  l3_isg_memory_cd4_t_cell               = "freq.L3_isg_memory_cd4_t_cell",
  l3_isg_memory_cd8_t_cell               = "freq.L3_isg_memory_cd8_t_cell",
  l3_isg_naive_b_cell                    = "freq.L3_isg_naive_b_cell",
  l3_isg_naive_cd4_t_cell                = "freq.L3_isg_naive_cd4_t_cell",
  l3_isg_naive_cd8_t_cell                = "freq.L3_isg_naive_cd8_t_cell",
  l3_klrf1_effector_vd1_gdt              = "freq.L3_klrf1_effector_vd1_gdt",
  l3_klrf1_effector_vd1_gdt_alt          = "freq.L3_klrf1_effector_vd1_gdt.1",
  l3_klrf1_gzmb_cd27_memory_cd4_t_cell   = "freq.L3_klrf1_gzmb_cd27_memory_cd4_t_cell",
  l3_naive_vd1_gdt                       = "freq.L3_naive_vd1_gdt",
  l3_pdc                                 = "freq.L3_pdc",
  l3_plasma_cell                         = "freq.L3_plasma_cell",
  l3_platelet                            = "freq.L3_platelet",
  l3_proliferating_nk_cell               = "freq.L3_proliferating_nk_cell",
  l3_proliferating_t_cell                = "freq.L3_proliferating_t_cell",
  l3_sox4_naive_cd4_t_cell               = "freq.L3_sox4_naive_cd4_t_cell",
  l3_sox4_naive_cd8_t_cell               = "freq.L3_sox4_naive_cd8_t_cell",
  l3_sox4_vd1_gdt                        = "freq.L3_sox4_vd1_gdt",
  l3_transitional_b_cell                 = "freq.L3_transitional_b_cell",
  l3_type_2_polarized_memory_b_cell      = "freq.L3_type_2_polarized_memory_b_cell",
  l3_adaptive_nk_cell                    = "freq.L3_adaptive_nk_cell",
  l3_clp_cell                            = "freq.L3_clp_cell",
  l3_cm_cd4_t_cell                       = "freq.L3_cm_cd4_t_cell",
  l3_cm_cd8_t_cell                       = "freq.L3_cm_cd8_t_cell",
  l3_cmp_cell                            = "freq.L3_cmp_cell",
  l3_gzmb_cd27_em_cd4_t_cell             = "freq.L3_gzmb_cd27_em_cd4_t_cell",
  l3_gzmb_cd27_em_cd4_t_cell_alt         = "freq.L3_gzmb_cd27_em_cd4_t_cell.1",
  l3_gzmk_cd27_em_cd8_t_cell             = "freq.L3_gzmk_cd27_em_cd8_t_cell",
  l3_gzmk_cd27_em_cd8_t_cell_alt         = "freq.L3_gzmk_cd27_em_cd8_t_cell.1",
  l3_gzmk_memory_cd4_treg                = "freq.L3_gzmk_memory_cd4_treg",
  l3_klrb1_memory_cd4_treg               = "freq.L3_klrb1_memory_cd4_treg",
  l3_klrb1_memory_cd8_treg               = "freq.L3_klrb1_memory_cd8_treg",
  l3_klrf1_gzmb_cd27_em_cd8_t_cell       = "freq.L3_klrf1_gzmb_cd27_em_cd8_t_cell",
  l3_klrf1_gzmb_cd27_em_cd8_t_cell_alt   = "freq.L3_klrf1_gzmb_cd27_em_cd8_t_cell.1",
  l3_memory_cd4_treg                     = "freq.L3_memory_cd4_treg",
  l3_memory_cd8_treg                     = "freq.L3_memory_cd8_treg",
  l3_naive_cd4_treg                      = "freq.L3_naive_cd4_treg",
  l3_baeomap_cell                        = "freq.L3_baeomap_cell"
)

# PARENT_MAP: child alias -> parent alias, used by 02_structure_learning_cv.Rmd
# to blacklist direct parent-child edges (both directions). Covers 99 of the
# 100 possible child nodes (29 L2 + 71 L3 - 1 unmapped baeomap_cell = 99).
#
# CONFIDENCE NOTE: entries below fall into two categories --
#   (a) name-substring matches (e.g. l2_cd14_monocyte contains "monocyte" ->
#       l1_monocyte) -- high confidence, mechanical
#   (b) standard immunology nomenclature not resolvable by substring alone
#       (e.g. "CM"/"EM" = central/effector memory, "CLP"/"CMP" = common
#       lymphoid/myeloid progenitor, Treg subsets, "adaptive NK") -- these
#       are domain-knowledge inferences and worth a second look from someone
#       with immunology expertise before treating them as ground truth:
#       l2_adaptive_nk_cell parent, l3_clp_cell, l3_cmp_cell, l3_cm_cd4_t_cell,
#       l3_cm_cd8_t_cell, l3_gzmb_cd27_em_cd4_t_cell(+_alt),
#       l3_gzmk_cd27_em_cd8_t_cell(+_alt), l3_klrf1_gzmb_cd27_em_cd8_t_cell(+_alt),
#       and all *_treg entries.
PARENT_MAP <- c(
  # L2 -> L1
  l2_asdc                 = "l1_dc",
  l2_cd14_monocyte        = "l1_monocyte",
  l2_cd16_monocyte        = "l1_monocyte",
  l2_cd56bright_nk_cell   = "l1_nk_cell",
  l2_cd56dim_nk_cell      = "l1_nk_cell",
  l2_cd8aa                = "l1_t_cell",
  l2_cdc1                 = "l1_dc",
  l2_cdc2                 = "l1_dc",
  l2_dn_t_cell            = "l1_t_cell",
  l2_effector_b_cell      = "l1_b_cell",
  l2_erythrocyte          = "l1_erythrocyte",
  l2_gdt                  = "l1_t_cell",
  l2_ilc                  = "l1_ilc",
  l2_intermediate_monocyte = "l1_monocyte",
  l2_mait                 = "l1_t_cell",
  l2_memory_b_cell        = "l1_b_cell",
  l2_memory_cd4_t_cell    = "l1_t_cell",
  l2_memory_cd8_t_cell    = "l1_t_cell",
  l2_naive_b_cell         = "l1_b_cell",
  l2_naive_cd4_t_cell     = "l1_t_cell",
  l2_naive_cd8_t_cell     = "l1_t_cell",
  l2_pdc                  = "l1_dc",
  l2_plasma_cell          = "l1_b_cell",
  l2_platelet             = "l1_platelet",
  l2_progenitor_cell      = "l1_progenitor_cell",
  l2_proliferating_nk_cell = "l1_nk_cell",
  l2_proliferating_t_cell = "l1_t_cell",
  l2_transitional_b_cell  = "l1_b_cell",
  l2_treg                 = "l1_t_cell",
  # L3 -> L2 (baeomap_cell omitted: no confident parent identified)
  l3_activated_memory_b_cell             = "l2_memory_b_cell",
  l3_asdc                                = "l2_asdc",
  l3_c1q_cd16_monocyte                   = "l2_cd16_monocyte",
  l3_cd14_cdc2                           = "l2_cdc2",
  l3_cd27_effector_b_cell                = "l2_effector_b_cell",
  l3_cd27_effector_b_cell_alt            = "l2_effector_b_cell",
  l3_cd4_mait                            = "l2_mait",
  l3_cd56bright_nk_cell                  = "l2_cd56bright_nk_cell",
  l3_cd8_mait                            = "l2_mait",
  l3_cd8aa                               = "l2_cd8aa",
  l3_cd95_memory_b_cell                  = "l2_memory_b_cell",
  l3_cdc1                                = "l2_cdc1",
  l3_core_cd14_monocyte                  = "l2_cd14_monocyte",
  l3_core_cd16_monocyte                  = "l2_cd16_monocyte",
  l3_core_memory_b_cell                  = "l2_memory_b_cell",
  l3_core_naive_b_cell                   = "l2_naive_b_cell",
  l3_core_naive_cd4_t_cell               = "l2_naive_cd4_t_cell",
  l3_core_naive_cd8_t_cell               = "l2_naive_cd8_t_cell",
  l3_dn_t_cell                           = "l2_dn_t_cell",
  l3_early_memory_b_cell                 = "l2_memory_b_cell",
  l3_erythrocyte                         = "l2_erythrocyte",
  l3_gzmb_vd2_gdt                        = "l2_gdt",
  l3_gzmk_cd56dim_nk_cell                = "l2_cd56dim_nk_cell",
  l3_gzmk_cd56dim_nk_cell_alt            = "l2_cd56dim_nk_cell",
  l3_gzmk_vd2_gdt                        = "l2_gdt",
  l3_hla_drhi_cdc2                       = "l2_cdc2",
  l3_il1b_cd14_monocyte                  = "l2_cd14_monocyte",
  l3_ilc                                 = "l2_ilc",
  l3_intermediate_monocyte               = "l2_intermediate_monocyte",
  l3_isg_cd14_monocyte                   = "l2_cd14_monocyte",
  l3_isg_cd16_monocyte                   = "l2_cd16_monocyte",
  l3_isg_cd56dim_nk_cell                 = "l2_cd56dim_nk_cell",
  l3_isg_cdc2                            = "l2_cdc2",
  l3_isg_mait                            = "l2_mait",
  l3_isg_memory_cd4_t_cell               = "l2_memory_cd4_t_cell",
  l3_isg_memory_cd8_t_cell               = "l2_memory_cd8_t_cell",
  l3_isg_naive_b_cell                    = "l2_naive_b_cell",
  l3_isg_naive_cd4_t_cell                = "l2_naive_cd4_t_cell",
  l3_isg_naive_cd8_t_cell                = "l2_naive_cd8_t_cell",
  l3_klrf1_effector_vd1_gdt              = "l2_gdt",
  l3_klrf1_effector_vd1_gdt_alt          = "l2_gdt",
  l3_klrf1_gzmb_cd27_memory_cd4_t_cell   = "l2_memory_cd4_t_cell",
  l3_naive_vd1_gdt                       = "l2_gdt",
  l3_pdc                                 = "l2_pdc",
  l3_plasma_cell                         = "l2_plasma_cell",
  l3_platelet                            = "l2_platelet",
  l3_proliferating_nk_cell               = "l2_proliferating_nk_cell",
  l3_proliferating_t_cell                = "l2_proliferating_t_cell",
  l3_sox4_naive_cd4_t_cell               = "l2_naive_cd4_t_cell",
  l3_sox4_naive_cd8_t_cell               = "l2_naive_cd8_t_cell",
  l3_sox4_vd1_gdt                        = "l2_gdt",
  l3_transitional_b_cell                 = "l2_transitional_b_cell",
  l3_type_2_polarized_memory_b_cell      = "l2_memory_b_cell",
  l3_adaptive_nk_cell                    = "l2_cd56dim_nk_cell",
  l3_clp_cell                            = "l2_progenitor_cell",
  l3_cm_cd4_t_cell                       = "l2_memory_cd4_t_cell",
  l3_cm_cd8_t_cell                       = "l2_memory_cd8_t_cell",
  l3_cmp_cell                            = "l2_progenitor_cell",
  l3_gzmb_cd27_em_cd4_t_cell             = "l2_memory_cd4_t_cell",
  l3_gzmb_cd27_em_cd4_t_cell_alt         = "l2_memory_cd4_t_cell",
  l3_gzmk_cd27_em_cd8_t_cell             = "l2_memory_cd8_t_cell",
  l3_gzmk_cd27_em_cd8_t_cell_alt         = "l2_memory_cd8_t_cell",
  l3_gzmk_memory_cd4_treg                = "l2_treg",
  l3_klrb1_memory_cd4_treg               = "l2_treg",
  l3_klrb1_memory_cd8_treg               = "l2_treg",
  l3_klrf1_gzmb_cd27_em_cd8_t_cell       = "l2_memory_cd8_t_cell",
  l3_klrf1_gzmb_cd27_em_cd8_t_cell_alt   = "l2_memory_cd8_t_cell",
  l3_memory_cd4_treg                     = "l2_treg",
  l3_memory_cd8_treg                     = "l2_treg",
  l3_naive_cd4_treg                      = "l2_treg"
)

NODE_MAP <- c(ROOT_MAP, CONTINUOUS_MAP)
ROOTS <- names(ROOT_MAP)

# Cell-frequency percentages -- no log transform applied here; revisit if
# diagnostics suggest right-skew (unlike clinical labs, these weren't
# pre-flagged as skewed in the original project handoff).
LOG_TRANSFORM <- character(0)

REFERENCE_LEVELS <- list(
  age_group = "Young Adult",
  sex       = "Female",
  cmv       = "Negative"
)
