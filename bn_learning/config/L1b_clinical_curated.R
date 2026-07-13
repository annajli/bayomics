# config/L1b_clinical_curated.R
# L1b: Clinical-only (curated) — 3 roots + 23 clinical continuous nodes
# Literature-curated immune-aging panel

NETWORK_ID <- "L1b_clinical_curated"
SEED <- 43L

EXCLUDE_CANDIDATES <- c(
  "clinical.chem.egfr_aa",
  "clinical.chem.egfr_non_aa",
  "clinical.infl.anti_ccp3",
  "clinical.infl.anti_ccp31",
  "clinical.infl.rf_iga_result",
  "clinical.infl.rf_igm_result"
)
DETERMINISTIC_EXCLUDE <- c("clinical.nlr")

ROOT_MAP <- c(
  age_group = "subject.ageGroup",
  sex       = "subject.biologicalSex",
  cmv       = "cmv.igg_serology_interpretation"
)

CONTINUOUS_MAP <- c(
  # Inflammation (2)
  crp              = "clinical.infl.hs_crp",
  esr              = "clinical.infl.esr",
  # Immune cells (4)
  lymphocytes      = "clinical.bc.lymphocyte_count",
  neutrophils      = "clinical.bc.neutrophil_count",
  monocytes        = "clinical.bc.monocyte_count",
  eosinophils      = "clinical.bc.eosinophil_count",
  # RBC / platelets (3)
  hemoglobin       = "clinical.bc.hemoglobin",
  platelets        = "clinical.bc.platelet_count",
  rdw              = "clinical.bc.rdw",
  # Lipids (3)
  hdl              = "clinical.lip.cholesterol_hdl",
  ldl              = "clinical.lip.cholesterol_ldl",
  triglycerides    = "clinical.lip.triglycerides",
  # Liver (3)
  alt              = "clinical.chem.alt",
  ast              = "clinical.chem.ast",
  albumin          = "clinical.chem.albumin",
  # Kidney (2)
  creatinine       = "clinical.chem.creatinine",
  bun              = "clinical.chem.bun",
  # Metabolic (2)
  glucose          = "clinical.chem.glucose",
  bmi              = "clinical.am.bmi",
  # CBC context (2)
  wbc              = "clinical.bc.wbc",
  hematocrit       = "clinical.bc.hematocrit",
  # Electrolytes (2)
  potassium        = "clinical.chem.potassium",
  calcium          = "clinical.chem.calcium"
)

NODE_MAP <- c(ROOT_MAP, CONTINUOUS_MAP)
ROOTS <- names(ROOT_MAP)

LOG_TRANSFORM <- c("crp", "esr", "alt", "ast", "triglycerides", "glucose")

REFERENCE_LEVELS <- list(
  age_group = "Young Adult",
  sex       = "Female",
  cmv       = "Negative"
)
