# config/L_all.R
# L_all: Kitchen-sink — ALL modalities combined
# 3 roots + all non-excluded continuous features from clinical, olink, freq,
# pseudobulk, and whole blood. Exact node count determined at data-prep time
# after missingness filtering and ALR transform on freq columns.
#
# EXPERIMENTAL: node-to-sample ratio will be ~8:1 (far beyond any validated
# network in this project — L6's 1.9 was the highest tested). maxp set
# aggressively low (5). Sponsor guidance: "not concerned with overfitting" —
# run and see what happens, but expect iteration. If catastrophic, try
# MAXP <- 3L. If stable, try MAXP <- 7L to check for missed structure.
#
# No whitelist used — let structure learning find edges on its own. Compare
# the learned structure against de_whitelist_candidates.csv post-hoc.
#
# CONTINUOUS_MAP is NOT defined here — it is built programmatically by
# 01_data_prep_all.Rmd from FEATURE_PREFIXES minus exclusions and
# missingness-filtered columns. NODE_MAP is also built at prep time.

NETWORK_ID <- "L_all"
SEED <- 50L
MAXP <- 5L

# --- Column discovery (used by 01_data_prep_all.Rmd) ---
FEATURE_PREFIXES <- c("clinical.", "olink.", "freq.", "pb.", "wb.")

# Per-column missingness threshold: columns with more than this fraction of
# samples missing are dropped before complete-case analysis. Set to 0.05
# (5%) to match L6's approach — with this many columns, even moderate per-
# column missingness compounds to unacceptable complete-case sample loss.
MISSINGNESS_THRESHOLD <- 0.05

# --- Exclusions (union of all individual-network exclusion lists) ---
EXCLUDE_CANDIDATES <- c(
  # Clinical (>20% missing)
  "clinical.chem.egfr_aa",
  "clinical.chem.egfr_non_aa",
  "clinical.infl.anti_ccp3",
  "clinical.infl.anti_ccp31",
  "clinical.infl.rf_iga_result",
  "clinical.infl.rf_igm_result",
  "clinical.am.bmi",
  "clinical.am.height",
  # Olink (>20% missing)
  "olink.MMP9",
  # WB QC flag (boolean, not biology)
  "wb.low_lib_flag"
)

DETERMINISTIC_EXCLUDE <- c(
  # Clinical deterministic/derived columns
  "clinical.nlr",
  "clinical.bc.mch",
  "clinical.bc.mchc",
  "clinical.bc.mcv",
  "clinical.bc.perc_basophils",
  "clinical.bc.perc_eosinophils",
  "clinical.bc.perc_lymphocytes",
  "clinical.bc.perc_monocytes",
  "clinical.bc.perc_neutrophils",
  "clinical.lip.cholesterol_non_hdl",
  "clinical.lip.chlesterol_hdl_ratio",
  # Freq exact duplicates (.2/.3 suffix variants)
  "freq.L3_cd27_effector_b_cell.2",              "freq.L3_cd27_effector_b_cell.3",
  "freq.L3_gzmb_cd27_em_cd4_t_cell.2",            "freq.L3_gzmb_cd27_em_cd4_t_cell.3",
  "freq.L3_gzmk_cd27_em_cd8_t_cell.2",            "freq.L3_gzmk_cd27_em_cd8_t_cell.3",
  "freq.L3_gzmk_cd56dim_nk_cell.2",               "freq.L3_gzmk_cd56dim_nk_cell.3",
  "freq.L3_klrf1_effector_vd1_gdt.2",             "freq.L3_klrf1_effector_vd1_gdt.3",
  "freq.L3_klrf1_gzmb_cd27_em_cd8_t_cell.2",      "freq.L3_klrf1_gzmb_cd27_em_cd8_t_cell.3"
)

ROOT_MAP <- c(
  age_group = "subject.ageGroup",
  sex       = "subject.biologicalSex",
  cmv       = "cmv.igg_serology_interpretation"
)
ROOTS <- names(ROOT_MAP)

# Right-skewed clinical variables that need log1p transform.
# Referenced by csv column names (not aliases) — the data prep notebook
# maps these to aliases after column discovery.
LOG_TRANSFORM_CSV <- c(
  "clinical.infl.hs_crp",
  "clinical.infl.esr",
  "clinical.chem.alt",
  "clinical.chem.ast",
  "clinical.lip.triglycerides",
  "clinical.chem.glucose"
)

REFERENCE_LEVELS <- list(
  age_group = "Young Adult",
  sex       = "Female",
  cmv       = "Negative"
)

# PARENT_MAP for freq ALR transform — identical to L5 config.
# Included here so 01_data_prep_all.Rmd can apply the ALR transform to
# the freq columns within the kitchen-sink network. Only the freq columns
# are affected; all other modalities pass through untransformed.
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
