# config/L6_pb_signaling.R
# Bootstrap pipeline version — node maps and exclusions identical to bn_cv.
# L6: PBMC functional/signaling response-only -- 3 roots + 144 nodes
# Phospho-flow cytometry: PBMCs stimulated with 6 different signals, measured
# for intracellular signaling pathway activation (not just cell counts --
# functional responsiveness) across 24 cell types. Original data had 65 cell
# types x 6 conditions (390 columns) + 65 QC flag columns (455 total); this
# config restricts to the 24 cell types complete enough to support real
# complete-case analysis (see missingness note below).

NETWORK_ID <- "L6_pb_signaling"
SEED <- 49L

# CELL-TYPE SELECTION: of 65 measured cell types, only the 24 below have
# <=5% per-subject missingness. The other 41 were dropped entirely (not
# just flagged) -- missingness here isn't random dropout, it's driven by
# cell-type rarity: if a subject didn't have enough cells of a given type,
# ALL 6 signaling readouts for that type are missing together for that
# subject (confirmed: 100% consistent within-cell-type missingness pattern
# across all 65 types). Combining many even-mildly-incomplete columns
# compounds fast -- at a 10% per-column cutoff, complete-case n crashes to
# 56/92; at 20%, to 18/92. The 5% cutoff keeps complete-case n at 81/92
# (12% sample loss), consistent with the <15% threshold already used
# elsewhere in this project to justify the bmi/height/MMP9 exclusions.
# Dropped cell types are biologically rare/small populations (ASDCs, common
# lymphoid progenitors, small Treg/MAIT subsets) -- excluded for the same
# data-availability reason as bmi/height, not a curation judgment call.
#
# QC FLAGS EXCLUDED: all 65 `*_low_ncells_flag` columns are boolean QC
# metadata (same category as wb.low_lib_flag in L4) -- not biology, not
# included in CONTINUOUS_MAP regardless of which cell types are kept.
EXCLUDE_CANDIDATES <- character(0)  # cell-type-level exclusion documented
                                      # above; no additional missingness
                                      # exclusions needed within the kept 24

# DETERMINISTIC_EXCLUDE: none. Checked pairwise correlations between the 6
# stimulation conditions within each of the 24 kept cell types (360 pairs
# total). 46 pairs exceed r=0.9 -- overwhelmingly ifn_alpha vs ifn_gamma
# (~20/24 cell types, r=0.90-0.98), consistent with both interferon types
# converging on shared STAT1 phosphorylation signaling, not an arithmetic
# identity. Max observed correlation was 0.976, well short of the near-1.0
# exact-duplicate pattern that justified exclusions in L1a/L5. Applying the
# same standard used for the clinical electrolyte cluster (real, moderate-
# to-high correlation from shared biology is left alone, not blacklisted,
# unless verified as a closed-form identity): no exclusions applied here.
# NOTE FOR INTERPRETATION: if ifn_alpha<->ifn_gamma edges appear frequently
# across many cell types in the resulting network, that reflects real shared
# pathway biology, not a modeling artifact.
DETERMINISTIC_EXCLUDE <- character(0)

ROOT_MAP <- c(
  age_group = "subject.ageGroup",
  sex       = "subject.biologicalSex",
  cmv       = "cmv.igg_serology_interpretation"
)

# Continuous PBMC signaling nodes (144 = 24 cell types x 6 stimulation conditions)
CONTINUOUS_MAP <- c(
  # --- cm_cd4_t_cell ---
  cm_cd4_t_cell_ifn_alpha                       = "pb.cm_cd4_t_cell_ifn_alpha",
  cm_cd4_t_cell_ifn_gamma                       = "pb.cm_cd4_t_cell_ifn_gamma",
  cm_cd4_t_cell_il2_stat5                       = "pb.cm_cd4_t_cell_il2_stat5",
  cm_cd4_t_cell_il6_jak_stat3                   = "pb.cm_cd4_t_cell_il6_jak_stat3",
  cm_cd4_t_cell_inflammatory                    = "pb.cm_cd4_t_cell_inflammatory",
  cm_cd4_t_cell_tnfa_nfkb                       = "pb.cm_cd4_t_cell_tnfa_nfkb",
  # --- cd56bright_nk_cell ---
  cd56bright_nk_cell_ifn_alpha                  = "pb.cd56bright_nk_cell_ifn_alpha",
  cd56bright_nk_cell_ifn_gamma                  = "pb.cd56bright_nk_cell_ifn_gamma",
  cd56bright_nk_cell_il2_stat5                  = "pb.cd56bright_nk_cell_il2_stat5",
  cd56bright_nk_cell_il6_jak_stat3              = "pb.cd56bright_nk_cell_il6_jak_stat3",
  cd56bright_nk_cell_inflammatory               = "pb.cd56bright_nk_cell_inflammatory",
  cd56bright_nk_cell_tnfa_nfkb                  = "pb.cd56bright_nk_cell_tnfa_nfkb",
  # --- gzmb_cd27_em_cd4_t_cell ---
  gzmb_cd27_em_cd4_t_cell_ifn_alpha             = "pb.gzmb_cd27_em_cd4_t_cell_ifn_alpha",
  gzmb_cd27_em_cd4_t_cell_ifn_gamma             = "pb.gzmb_cd27_em_cd4_t_cell_ifn_gamma",
  gzmb_cd27_em_cd4_t_cell_il2_stat5             = "pb.gzmb_cd27_em_cd4_t_cell_il2_stat5",
  gzmb_cd27_em_cd4_t_cell_il6_jak_stat3         = "pb.gzmb_cd27_em_cd4_t_cell_il6_jak_stat3",
  gzmb_cd27_em_cd4_t_cell_inflammatory          = "pb.gzmb_cd27_em_cd4_t_cell_inflammatory",
  gzmb_cd27_em_cd4_t_cell_tnfa_nfkb             = "pb.gzmb_cd27_em_cd4_t_cell_tnfa_nfkb",
  # --- core_naive_cd8_t_cell ---
  core_naive_cd8_t_cell_ifn_alpha               = "pb.core_naive_cd8_t_cell_ifn_alpha",
  core_naive_cd8_t_cell_ifn_gamma               = "pb.core_naive_cd8_t_cell_ifn_gamma",
  core_naive_cd8_t_cell_il2_stat5               = "pb.core_naive_cd8_t_cell_il2_stat5",
  core_naive_cd8_t_cell_il6_jak_stat3           = "pb.core_naive_cd8_t_cell_il6_jak_stat3",
  core_naive_cd8_t_cell_inflammatory            = "pb.core_naive_cd8_t_cell_inflammatory",
  core_naive_cd8_t_cell_tnfa_nfkb               = "pb.core_naive_cd8_t_cell_tnfa_nfkb",
  # --- core_naive_b_cell ---
  core_naive_b_cell_ifn_alpha                   = "pb.core_naive_b_cell_ifn_alpha",
  core_naive_b_cell_ifn_gamma                   = "pb.core_naive_b_cell_ifn_gamma",
  core_naive_b_cell_il2_stat5                   = "pb.core_naive_b_cell_il2_stat5",
  core_naive_b_cell_il6_jak_stat3               = "pb.core_naive_b_cell_il6_jak_stat3",
  core_naive_b_cell_inflammatory                = "pb.core_naive_b_cell_inflammatory",
  core_naive_b_cell_tnfa_nfkb                   = "pb.core_naive_b_cell_tnfa_nfkb",
  # --- core_naive_cd4_t_cell ---
  core_naive_cd4_t_cell_ifn_alpha               = "pb.core_naive_cd4_t_cell_ifn_alpha",
  core_naive_cd4_t_cell_ifn_gamma               = "pb.core_naive_cd4_t_cell_ifn_gamma",
  core_naive_cd4_t_cell_il2_stat5               = "pb.core_naive_cd4_t_cell_il2_stat5",
  core_naive_cd4_t_cell_il6_jak_stat3           = "pb.core_naive_cd4_t_cell_il6_jak_stat3",
  core_naive_cd4_t_cell_inflammatory            = "pb.core_naive_cd4_t_cell_inflammatory",
  core_naive_cd4_t_cell_tnfa_nfkb               = "pb.core_naive_cd4_t_cell_tnfa_nfkb",
  # --- core_memory_b_cell ---
  core_memory_b_cell_ifn_alpha                  = "pb.core_memory_b_cell_ifn_alpha",
  core_memory_b_cell_ifn_gamma                  = "pb.core_memory_b_cell_ifn_gamma",
  core_memory_b_cell_il2_stat5                  = "pb.core_memory_b_cell_il2_stat5",
  core_memory_b_cell_il6_jak_stat3              = "pb.core_memory_b_cell_il6_jak_stat3",
  core_memory_b_cell_inflammatory               = "pb.core_memory_b_cell_inflammatory",
  core_memory_b_cell_tnfa_nfkb                  = "pb.core_memory_b_cell_tnfa_nfkb",
  # --- core_cd14_monocyte ---
  core_cd14_monocyte_ifn_alpha                  = "pb.core_cd14_monocyte_ifn_alpha",
  core_cd14_monocyte_ifn_gamma                  = "pb.core_cd14_monocyte_ifn_gamma",
  core_cd14_monocyte_il2_stat5                  = "pb.core_cd14_monocyte_il2_stat5",
  core_cd14_monocyte_il6_jak_stat3              = "pb.core_cd14_monocyte_il6_jak_stat3",
  core_cd14_monocyte_inflammatory               = "pb.core_cd14_monocyte_inflammatory",
  core_cd14_monocyte_tnfa_nfkb                  = "pb.core_cd14_monocyte_tnfa_nfkb",
  # --- gzmk_cd56dim_nk_cell ---
  gzmk_cd56dim_nk_cell_ifn_alpha                = "pb.gzmk_cd56dim_nk_cell_ifn_alpha",
  gzmk_cd56dim_nk_cell_ifn_gamma                = "pb.gzmk_cd56dim_nk_cell_ifn_gamma",
  gzmk_cd56dim_nk_cell_il2_stat5                = "pb.gzmk_cd56dim_nk_cell_il2_stat5",
  gzmk_cd56dim_nk_cell_il6_jak_stat3            = "pb.gzmk_cd56dim_nk_cell_il6_jak_stat3",
  gzmk_cd56dim_nk_cell_inflammatory             = "pb.gzmk_cd56dim_nk_cell_inflammatory",
  gzmk_cd56dim_nk_cell_tnfa_nfkb                = "pb.gzmk_cd56dim_nk_cell_tnfa_nfkb",
  # --- naive_cd4_treg ---
  naive_cd4_treg_ifn_alpha                      = "pb.naive_cd4_treg_ifn_alpha",
  naive_cd4_treg_ifn_gamma                      = "pb.naive_cd4_treg_ifn_gamma",
  naive_cd4_treg_il2_stat5                      = "pb.naive_cd4_treg_il2_stat5",
  naive_cd4_treg_il6_jak_stat3                  = "pb.naive_cd4_treg_il6_jak_stat3",
  naive_cd4_treg_inflammatory                   = "pb.naive_cd4_treg_inflammatory",
  naive_cd4_treg_tnfa_nfkb                      = "pb.naive_cd4_treg_tnfa_nfkb",
  # --- transitional_b_cell ---
  transitional_b_cell_ifn_alpha                 = "pb.transitional_b_cell_ifn_alpha",
  transitional_b_cell_ifn_gamma                 = "pb.transitional_b_cell_ifn_gamma",
  transitional_b_cell_il2_stat5                 = "pb.transitional_b_cell_il2_stat5",
  transitional_b_cell_il6_jak_stat3             = "pb.transitional_b_cell_il6_jak_stat3",
  transitional_b_cell_inflammatory              = "pb.transitional_b_cell_inflammatory",
  transitional_b_cell_tnfa_nfkb                 = "pb.transitional_b_cell_tnfa_nfkb",
  # --- isg_cd14_monocyte ---
  isg_cd14_monocyte_ifn_alpha                   = "pb.isg_cd14_monocyte_ifn_alpha",
  isg_cd14_monocyte_ifn_gamma                   = "pb.isg_cd14_monocyte_ifn_gamma",
  isg_cd14_monocyte_il2_stat5                   = "pb.isg_cd14_monocyte_il2_stat5",
  isg_cd14_monocyte_il6_jak_stat3               = "pb.isg_cd14_monocyte_il6_jak_stat3",
  isg_cd14_monocyte_inflammatory                = "pb.isg_cd14_monocyte_inflammatory",
  isg_cd14_monocyte_tnfa_nfkb                   = "pb.isg_cd14_monocyte_tnfa_nfkb",
  # --- memory_cd4_treg ---
  memory_cd4_treg_ifn_alpha                     = "pb.memory_cd4_treg_ifn_alpha",
  memory_cd4_treg_ifn_gamma                     = "pb.memory_cd4_treg_ifn_gamma",
  memory_cd4_treg_il2_stat5                     = "pb.memory_cd4_treg_il2_stat5",
  memory_cd4_treg_il6_jak_stat3                 = "pb.memory_cd4_treg_il6_jak_stat3",
  memory_cd4_treg_inflammatory                  = "pb.memory_cd4_treg_inflammatory",
  memory_cd4_treg_tnfa_nfkb                     = "pb.memory_cd4_treg_tnfa_nfkb",
  # --- core_cd16_monocyte ---
  core_cd16_monocyte_ifn_alpha                  = "pb.core_cd16_monocyte_ifn_alpha",
  core_cd16_monocyte_ifn_gamma                  = "pb.core_cd16_monocyte_ifn_gamma",
  core_cd16_monocyte_il2_stat5                  = "pb.core_cd16_monocyte_il2_stat5",
  core_cd16_monocyte_il6_jak_stat3              = "pb.core_cd16_monocyte_il6_jak_stat3",
  core_cd16_monocyte_inflammatory               = "pb.core_cd16_monocyte_inflammatory",
  core_cd16_monocyte_tnfa_nfkb                  = "pb.core_cd16_monocyte_tnfa_nfkb",
  # --- cd14_cdc2 ---
  cd14_cdc2_ifn_alpha                           = "pb.cd14_cdc2_ifn_alpha",
  cd14_cdc2_ifn_gamma                           = "pb.cd14_cdc2_ifn_gamma",
  cd14_cdc2_il2_stat5                           = "pb.cd14_cdc2_il2_stat5",
  cd14_cdc2_il6_jak_stat3                       = "pb.cd14_cdc2_il6_jak_stat3",
  cd14_cdc2_inflammatory                        = "pb.cd14_cdc2_inflammatory",
  cd14_cdc2_tnfa_nfkb                           = "pb.cd14_cdc2_tnfa_nfkb",
  # --- cm_cd8_t_cell ---
  cm_cd8_t_cell_ifn_alpha                       = "pb.cm_cd8_t_cell_ifn_alpha",
  cm_cd8_t_cell_ifn_gamma                       = "pb.cm_cd8_t_cell_ifn_gamma",
  cm_cd8_t_cell_il2_stat5                       = "pb.cm_cd8_t_cell_il2_stat5",
  cm_cd8_t_cell_il6_jak_stat3                   = "pb.cm_cd8_t_cell_il6_jak_stat3",
  cm_cd8_t_cell_inflammatory                    = "pb.cm_cd8_t_cell_inflammatory",
  cm_cd8_t_cell_tnfa_nfkb                       = "pb.cm_cd8_t_cell_tnfa_nfkb",
  # --- cd8_mait ---
  cd8_mait_ifn_alpha                            = "pb.cd8_mait_ifn_alpha",
  cd8_mait_ifn_gamma                            = "pb.cd8_mait_ifn_gamma",
  cd8_mait_il2_stat5                            = "pb.cd8_mait_il2_stat5",
  cd8_mait_il6_jak_stat3                        = "pb.cd8_mait_il6_jak_stat3",
  cd8_mait_inflammatory                         = "pb.cd8_mait_inflammatory",
  cd8_mait_tnfa_nfkb                            = "pb.cd8_mait_tnfa_nfkb",
  # --- intermediate_monocyte ---
  intermediate_monocyte_ifn_alpha               = "pb.intermediate_monocyte_ifn_alpha",
  intermediate_monocyte_ifn_gamma               = "pb.intermediate_monocyte_ifn_gamma",
  intermediate_monocyte_il2_stat5               = "pb.intermediate_monocyte_il2_stat5",
  intermediate_monocyte_il6_jak_stat3           = "pb.intermediate_monocyte_il6_jak_stat3",
  intermediate_monocyte_inflammatory            = "pb.intermediate_monocyte_inflammatory",
  intermediate_monocyte_tnfa_nfkb               = "pb.intermediate_monocyte_tnfa_nfkb",
  # --- hla_drhi_cdc2 ---
  hla_drhi_cdc2_ifn_alpha                       = "pb.hla_drhi_cdc2_ifn_alpha",
  hla_drhi_cdc2_ifn_gamma                       = "pb.hla_drhi_cdc2_ifn_gamma",
  hla_drhi_cdc2_il2_stat5                       = "pb.hla_drhi_cdc2_il2_stat5",
  hla_drhi_cdc2_il6_jak_stat3                   = "pb.hla_drhi_cdc2_il6_jak_stat3",
  hla_drhi_cdc2_inflammatory                    = "pb.hla_drhi_cdc2_inflammatory",
  hla_drhi_cdc2_tnfa_nfkb                       = "pb.hla_drhi_cdc2_tnfa_nfkb",
  # --- pdc ---
  pdc_ifn_alpha                                 = "pb.pdc_ifn_alpha",
  pdc_ifn_gamma                                 = "pb.pdc_ifn_gamma",
  pdc_il2_stat5                                 = "pb.pdc_il2_stat5",
  pdc_il6_jak_stat3                             = "pb.pdc_il6_jak_stat3",
  pdc_inflammatory                              = "pb.pdc_inflammatory",
  pdc_tnfa_nfkb                                 = "pb.pdc_tnfa_nfkb",
  # --- klrb1_memory_cd4_treg ---
  klrb1_memory_cd4_treg_ifn_alpha               = "pb.klrb1_memory_cd4_treg_ifn_alpha",
  klrb1_memory_cd4_treg_ifn_gamma               = "pb.klrb1_memory_cd4_treg_ifn_gamma",
  klrb1_memory_cd4_treg_il2_stat5               = "pb.klrb1_memory_cd4_treg_il2_stat5",
  klrb1_memory_cd4_treg_il6_jak_stat3           = "pb.klrb1_memory_cd4_treg_il6_jak_stat3",
  klrb1_memory_cd4_treg_inflammatory            = "pb.klrb1_memory_cd4_treg_inflammatory",
  klrb1_memory_cd4_treg_tnfa_nfkb               = "pb.klrb1_memory_cd4_treg_tnfa_nfkb",
  # --- gzmk_vd2_gdt ---
  gzmk_vd2_gdt_ifn_alpha                        = "pb.gzmk_vd2_gdt_ifn_alpha",
  gzmk_vd2_gdt_ifn_gamma                        = "pb.gzmk_vd2_gdt_ifn_gamma",
  gzmk_vd2_gdt_il2_stat5                        = "pb.gzmk_vd2_gdt_il2_stat5",
  gzmk_vd2_gdt_il6_jak_stat3                    = "pb.gzmk_vd2_gdt_il6_jak_stat3",
  gzmk_vd2_gdt_inflammatory                     = "pb.gzmk_vd2_gdt_inflammatory",
  gzmk_vd2_gdt_tnfa_nfkb                        = "pb.gzmk_vd2_gdt_tnfa_nfkb",
  # --- klrf1_gzmb_cd27_em_cd8_t_cell ---
  klrf1_gzmb_cd27_em_cd8_t_cell_ifn_alpha       = "pb.klrf1_gzmb_cd27_em_cd8_t_cell_ifn_alpha",
  klrf1_gzmb_cd27_em_cd8_t_cell_ifn_gamma       = "pb.klrf1_gzmb_cd27_em_cd8_t_cell_ifn_gamma",
  klrf1_gzmb_cd27_em_cd8_t_cell_il2_stat5       = "pb.klrf1_gzmb_cd27_em_cd8_t_cell_il2_stat5",
  klrf1_gzmb_cd27_em_cd8_t_cell_il6_jak_stat3   = "pb.klrf1_gzmb_cd27_em_cd8_t_cell_il6_jak_stat3",
  klrf1_gzmb_cd27_em_cd8_t_cell_inflammatory    = "pb.klrf1_gzmb_cd27_em_cd8_t_cell_inflammatory",
  klrf1_gzmb_cd27_em_cd8_t_cell_tnfa_nfkb       = "pb.klrf1_gzmb_cd27_em_cd8_t_cell_tnfa_nfkb",
  # --- isg_memory_cd4_t_cell ---
  isg_memory_cd4_t_cell_ifn_alpha               = "pb.isg_memory_cd4_t_cell_ifn_alpha",
  isg_memory_cd4_t_cell_ifn_gamma               = "pb.isg_memory_cd4_t_cell_ifn_gamma",
  isg_memory_cd4_t_cell_il2_stat5               = "pb.isg_memory_cd4_t_cell_il2_stat5",
  isg_memory_cd4_t_cell_il6_jak_stat3           = "pb.isg_memory_cd4_t_cell_il6_jak_stat3",
  isg_memory_cd4_t_cell_inflammatory            = "pb.isg_memory_cd4_t_cell_inflammatory",
  isg_memory_cd4_t_cell_tnfa_nfkb               = "pb.isg_memory_cd4_t_cell_tnfa_nfkb"
)

NODE_MAP <- c(ROOT_MAP, CONTINUOUS_MAP)
ROOTS <- names(ROOT_MAP)

# Signaling readouts are already normalized/standardized -- checked mean
# ~=0 and presence of negative values (e.g. pb.cm_cd4_t_cell_ifn_alpha:
# mean=-3e-8, range -1.05 to 1.39), confirming this is not raw untransformed
# data. Could NOT confirm which specific transform was applied (arcsinh,
# z-score, or something else) from the CSV alone -- verify against the data
# source before assuming. No additional log transform applied here, since
# adding one on top of an already-normalized (and partly negative-valued)
# variable would be inappropriate regardless of which transform was used.
LOG_TRANSFORM <- character(0)

REFERENCE_LEVELS <- list(
  age_group = "Young Adult",
  sex       = "Female",
  cmv       = "Negative"
)
