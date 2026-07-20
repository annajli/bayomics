# config/L1a_clinical_full_v3.R
# Bootstrap pipeline version.
# L1a v3: Clinical-only (full) — 3 roots + 34 clinical continuous nodes
# Based on v2 (which excluded BMI/height for missingness), now also excluding
# 10 deterministic/derived columns that are arithmetic functions of other
# included nodes. These produced trivial edges (e.g. basophils → perc_basophils,
# hemoglobin ↔ hematocrit ↔ mch/mchc/mcv) that dominated the DAG and masked
# biologically interesting structure.

NETWORK_ID <- "L1a_clinical_full_v3"
SEED <- 42L

# Columns to drop (>20% missing)
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

# Deterministic / derived columns — arithmetic functions of other included
# nodes, producing trivial edges that dominate the DAG:
#   NLR = neutrophil_count / lymphocyte_count
#   mch = hemoglobin / rbc, mchc = hemoglobin / hematocrit, mcv = hematocrit / rbc
#   perc_X = X_count / wbc * 100 (5 columns)
#   non_hdl = total_chol - hdl, hdl_ratio = total_chol / hdl
DETERMINISTIC_EXCLUDE <- c(
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
  "clinical.lip.chlesterol_hdl_ratio"
)

# Short alias -> bn_ready_baseline.csv column name
# Roots (discrete)
ROOT_MAP <- c(
  age_group = "subject.ageGroup",
  sex       = "subject.biologicalSex",
  cmv       = "cmv.igg_serology_interpretation"
)

# Continuous clinical nodes (34 = 53 - 8 exclude_candidates - 11 deterministic)
CONTINUOUS_MAP <- c(
  weight           = "clinical.am.weight",
  basophils        = "clinical.bc.basophil_count",
  eosinophils      = "clinical.bc.eosinophil_count",
  hematocrit       = "clinical.bc.hematocrit",
  hemoglobin       = "clinical.bc.hemoglobin",
  lymphocytes      = "clinical.bc.lymphocyte_count",
  monocytes        = "clinical.bc.monocyte_count",
  mpv              = "clinical.bc.mpv",
  neutrophils      = "clinical.bc.neutrophil_count",
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
  hdl              = "clinical.lip.cholesterol_hdl",
  ldl              = "clinical.lip.cholesterol_ldl",
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
