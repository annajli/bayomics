# config/L4_wb_pathways.R
# Bootstrap pipeline version — node maps and exclusions identical to bn_cv.
# L4: Whole blood-only — 3 roots + 50 MSigDB Hallmark pathway scores
# Bulk RNA-seq from whole blood, scored against the 50 Hallmark gene sets.
# Each node is a single per-subject summary score (not individual gene
# expression) — a coarse, interpretable readout of pathway-level activity.

NETWORK_ID <- "L4_wb_pathways"
SEED <- 45L

# wb.low_lib_flag EXCLUDED — this is a boolean RNA-seq QC flag (low library
# size: 90 False / 1 True / 1 NA), not a Hallmark pathway score. It doesn't
# match the naming pattern of the other 50 columns and is categorical
# metadata about sequencing quality, not biology. Including it in
# CONTINUOUS_MAP would be a type error, not just noise.
EXCLUDE_CANDIDATES <- c("wb.low_lib_flag")
DETERMINISTIC_EXCLUDE <- character(0)  # Hallmark gene sets can share member
# genes and correlate biologically, but are not exact arithmetic functions
# of one another the way clinical.nlr or bc.mchc are — no exclusions here
# unless verified otherwise against real data.

# Missingness: all 50 real pathway columns are missing for the same 2/92
# subjects (2.2% each, uniformly) — consistent with 2 subjects having a
# single failed/missing RNA-seq run, not scattered per-column dropout.
# Complete-case n = 90/92. No individual-column exclusions needed.

ROOT_MAP <- c(
  age_group = "subject.ageGroup",
  sex       = "subject.biologicalSex",
  cmv       = "cmv.igg_serology_interpretation"
)

# Continuous whole-blood pathway nodes (50)
CONTINUOUS_MAP <- c(
  # Proliferation / cell cycle
  e2f_targets                          = "wb.e2f_targets",
  g2_m_checkpoint                      = "wb.g2_m_checkpoint",
  mitotic_spindle                      = "wb.mitotic_spindle",
  myc_targets_v1                       = "wb.myc_targets_v1",
  myc_targets_v2                       = "wb.myc_targets_v2",
  # DNA damage / stress response
  dna_repair                           = "wb.dna_repair",
  p53_pathway                          = "wb.p53_pathway",
  uv_response_dn                       = "wb.uv_response_dn",
  uv_response_up                       = "wb.uv_response_up",
  unfolded_protein_response            = "wb.unfolded_protein_response",
  reactive_oxygen_species_pathway      = "wb.reactive_oxygen_species_pathway",
  hypoxia                              = "wb.hypoxia",
  # Immune / inflammatory signaling
  inflammatory_response                = "wb.inflammatory_response",
  interferon_alpha_response            = "wb.interferon_alpha_response",
  interferon_gamma_response            = "wb.interferon_gamma_response",
  tnf_alpha_signaling_via_nf_kb        = "wb.tnf_alpha_signaling_via_nf_kb",
  il_2_stat5_signaling                 = "wb.il_2_stat5_signaling",
  il_6_jak_stat3_signaling             = "wb.il_6_jak_stat3_signaling",
  complement                           = "wb.complement",
  coagulation                          = "wb.coagulation",
  allograft_rejection                  = "wb.allograft_rejection",
  # Signaling pathways
  pi3k_akt_mtor_signaling              = "wb.pi3k_akt_mtor_signaling",
  mtorc1_signaling                     = "wb.mtorc1_signaling",
  kras_signaling_up                    = "wb.kras_signaling_up",
  kras_signaling_dn                    = "wb.kras_signaling_dn",
  tgf_beta_signaling                   = "wb.tgf_beta_signaling",
  wnt_beta_catenin_signaling           = "wb.wnt_beta_catenin_signaling",
  notch_signaling                      = "wb.notch_signaling",
  hedgehog_signaling                   = "wb.hedgehog_signaling",
  # Metabolic
  glycolysis                           = "wb.glycolysis",
  oxidative_phosphorylation            = "wb.oxidative_phosphorylation",
  fatty_acid_metabolism                = "wb.fatty_acid_metabolism",
  bile_acid_metabolism                 = "wb.bile_acid_metabolism",
  cholesterol_homeostasis              = "wb.cholesterol_homeostasis",
  heme_metabolism                      = "wb.heme_metabolism",
  xenobiotic_metabolism                = "wb.xenobiotic_metabolism",
  pperoxisome                          = "wb.pperoxisome",
  # Development / differentiation
  adipogenesis                         = "wb.adipogenesis",
  myogenesis                           = "wb.myogenesis",
  spermatogenesis                      = "wb.spermatogenesis",
  pancreas_beta_cells                  = "wb.pancreas_beta_cells",
  epithelial_mesenchymal_transition    = "wb.epithelial_mesenchymal_transition",
  angiogenesis                         = "wb.angiogenesis",
  # Hormone response
  androgen_response                    = "wb.androgen_response",
  estrogen_response_early              = "wb.estrogen_response_early",
  estrogen_response_late               = "wb.estrogen_response_late",
  # Cellular component / structure
  apical_junction                      = "wb.apical_junction",
  apical_surface                       = "wb.apical_surface",
  # Other
  apoptosis                            = "wb.apoptosis",
  protein_secretion                    = "wb.protein_secretion"
)

NODE_MAP <- c(ROOT_MAP, CONTINUOUS_MAP)
ROOTS <- names(ROOT_MAP)

# Hallmark scores are already normalized enrichment scores — no additional
# log transform needed (same reasoning as Olink's NPX values in L2).
LOG_TRANSFORM <- character(0)

REFERENCE_LEVELS <- list(
  age_group = "Young Adult",
  sex       = "Female",
  cmv       = "Negative"
)
