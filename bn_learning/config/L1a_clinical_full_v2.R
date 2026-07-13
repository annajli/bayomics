# config/L1a_clinical_full.R
# L1a: Clinical-only (full) — 3 roots + 47 clinical continuous nodes
# 53 clinical features minus 6 exclude_candidates (>20% missing)

NETWORK_ID <- "L1a_clinical_full"
SEED <- 42L

# Columns to drop (>20% missing)
# bmi/height added 2026-07: each 18/92 (19.6%) missing; complete-case analysis
# requiring all 47 original nodes dropped 23/92 (25%) samples, exceeding the
# 15% threshold. Consistent with L1b, which already excludes bmi for the same
# reason.
EXCLUDE_CANDIDATES <- c(
  "clinical.chem.egfr_aa",
  "clinical.chem.egfr_non_aa",
  "clinical.infl.anti_ccp3",
  "clinical.infl.anti_ccp31",
  "clinical.infl.rf_iga_result",
  "clinical.infl.rf_igm_result",
  "clinical.am.bmi",
  "clinical.am.height"
)

# Also exclude NLR (deterministic function of neutrophil/lymphocyte counts)
DETERMINISTIC_EXCLUDE <- c("clinical.nlr")

# Short alias -> bn_ready_baseline.csv column name
# Roots (discrete)
ROOT_MAP <- c(
  age_group = "subject.ageGroup",
  sex       = "subject.biologicalSex",
  cmv       = "cmv.igg_serology_interpretation"
)

# Continuous clinical nodes (47 = 53 - 6 exclude_candidates)
CONTINUOUS_MAP <- c(
  bmi              = "clinical.am.bmi",
  height           = "clinical.am.height",
  weight           = "clinical.am.weight",
  basophils        = "clinical.bc.basophil_count",
  eosinophils      = "clinical.bc.eosinophil_count",
  hematocrit       = "clinical.bc.hematocrit",
  hemoglobin       = "clinical.bc.hemoglobin",
  lymphocytes      = "clinical.bc.lymphocyte_count",
  mch              = "clinical.bc.mch",
  mchc             = "clinical.bc.mchc",
  mcv              = "clinical.bc.mcv",
  monocytes        = "clinical.bc.monocyte_count",
  mpv              = "clinical.bc.mpv",
  neutrophils      = "clinical.bc.neutrophil_count",
  perc_basophils   = "clinical.bc.perc_basophils",
  perc_eosinophils = "clinical.bc.perc_eosinophils",
  perc_lymphocytes = "clinical.bc.perc_lymphocytes",
  perc_monocytes   = "clinical.bc.perc_monocytes",
  perc_neutrophils = "clinical.bc.perc_neutrophils",
  platelets        = "clinical.bc.platelet_count",
  rdw              = "clinical.bc.rdw",
  rbc              = "clinical.bc.red_blood_cell_count",
  wbc              = "clinical.bc.wbc",
  albumin          = "clinical.chem.albumin",
  alk_phos         = "clinical.chem.alkaline_phosphatase",
  alt              = "clinical.chem.alt",
  ast              = "clinical.chem.ast",
  bun              = "clinical.chem.bun",
  calcium          = "clinical.chem.calcium",
  chloride         = "clinical.chem.cl",
  co2              = "clinical.chem.co2",
  creatinine       = "clinical.chem.creatinine",
  globulin         = "clinical.chem.globin",
  glucose          = "clinical.chem.glucose",
  potassium        = "clinical.chem.potassium",
  protein          = "clinical.chem.protein",
  sodium           = "clinical.chem.sodium",
  t_bili           = "clinical.chem.t_bili",
  esr              = "clinical.infl.esr",
  crp              = "clinical.infl.hs_crp",
  hdl_ratio        = "clinical.lip.chlesterol_hdl_ratio",
  hdl              = "clinical.lip.cholesterol_hdl",
  ldl              = "clinical.lip.cholesterol_ldl",
  non_hdl          = "clinical.lip.cholesterol_non_hdl",
  total_chol       = "clinical.lip.cholesterol_total",
  triglycerides    = "clinical.lip.triglycerides"
)

NODE_MAP <- c(ROOT_MAP, CONTINUOUS_MAP)
ROOTS <- names(ROOT_MAP)

# Variables to log1p-transform (right-skewed)
LOG_TRANSFORM <- c("crp", "esr", "alt", "ast", "triglycerides", "glucose")

# Factor reference levels for roots
REFERENCE_LEVELS <- list(
  age_group = "Young Adult",
  sex       = "Female",
  cmv       = "Negative"
)
