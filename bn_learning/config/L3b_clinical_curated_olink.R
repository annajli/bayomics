# config/L3b_clinical_curated_olink.R
# Bootstrap pipeline version — node maps and exclusions identical to bn_cv.
# L3b: Combined — curated (literature-based) clinical panel + Olink proteomics
# 3 roots + 22 clinical continuous + 34 Olink continuous = 59 nodes total
#
# Built by merging the current L1b_clinical_curated_v2.R and L2_olink_v2.R
# configs. Both source configs' exclusion decisions are preserved as-is:
#   - Clinical: 8 EXCLUDE_CANDIDATES (same missingness list as L1a/L1b) +
#     1 DETERMINISTIC_EXCLUDE (NLR only -- L1b's curated panel was never
#     reviewed for the broader deterministic-column issue found in L1a;
#     see Open Items in the project README)
#   - Olink: 1 EXCLUDE_CANDIDATE (MMP9, missingness)
# No new exclusions introduced by combining -- verified no alias name
# collisions between the 22 clinical and 34 Olink node names.
#
# NOT YET CHECKED: cross-modality collinearity between the curated clinical
# panel and Olink, same caveat as L3a.

NETWORK_ID <- "L3b_clinical_curated_olink"
SEED <- 48L

EXCLUDE_CANDIDATES <- c(
  # -- from L1b_clinical_curated_v2.R --
  "clinical.chem.egfr_aa",
  "clinical.chem.egfr_non_aa",
  "clinical.infl.anti_ccp3",
  "clinical.infl.anti_ccp31",
  "clinical.infl.rf_iga_result",
  "clinical.infl.rf_igm_result",
  "clinical.am.bmi",
  "clinical.am.height",
  # -- from L2_olink_v2.R --
  "olink.MMP9"
)

DETERMINISTIC_EXCLUDE <- c("clinical.nlr")  # from L1b_clinical_curated_v2.R only

ROOT_MAP <- c(
  age_group = "subject.ageGroup",
  sex       = "subject.biologicalSex",
  cmv       = "cmv.igg_serology_interpretation"
)

# Continuous nodes (56 = 22 clinical + 34 Olink)
CONTINUOUS_MAP <- c(
  # === Clinical (22) === (from L1b_clinical_curated_v2.R)
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
  # Metabolic (1)
  glucose          = "clinical.chem.glucose",
  # CBC context (2)
  wbc              = "clinical.bc.wbc",
  hematocrit       = "clinical.bc.hematocrit",
  # Electrolytes (2)
  potassium        = "clinical.chem.potassium",
  calcium          = "clinical.chem.calcium",
  # === Olink (34) === (from L2_olink_v2.R)
  IL6     = "olink.IL6",
  TNF     = "olink.TNF",
  CXCL8   = "olink.CXCL8",
  CXCL9   = "olink.CXCL9",
  CXCL10  = "olink.CXCL10",
  CXCL11  = "olink.CXCL11",
  IFNG    = "olink.IFNG",
  IL18    = "olink.IL18",
  IL10    = "olink.IL10",
  IL1RN   = "olink.IL1RN",
  IL7     = "olink.IL7",
  IL15    = "olink.IL15",
  CCL2    = "olink.CCL2",
  CCL3    = "olink.CCL3",
  CCL4    = "olink.CCL4",
  CCL8    = "olink.CCL8",
  CCL11   = "olink.CCL11",
  GDF15   = "olink.GDF15",
  FGF21   = "olink.FGF21",
  CHI3L1  = "olink.CHI3L1",
  IGFBP7  = "olink.IGFBP7",
  TNFRSF1A = "olink.TNFRSF1A",
  TNFRSF1B = "olink.TNFRSF1B",
  CD8A    = "olink.CD8A",
  GZMA    = "olink.GZMA",
  GZMB    = "olink.GZMB",
  KLRD1   = "olink.KLRD1",
  MMP12   = "olink.MMP12",
  VEGFA   = "olink.VEGFA",
  PDCD1   = "olink.PDCD1",
  LAG3    = "olink.LAG3",
  HAVCR2  = "olink.HAVCR2",
  CD274   = "olink.CD274",
  LEP     = "olink.LEP"
)

NODE_MAP <- c(ROOT_MAP, CONTINUOUS_MAP)
ROOTS <- names(ROOT_MAP)

# Log transform: clinical's right-skewed variables only (crp, esr, alt, ast,
# triglycerides all present in L1b's curated panel; glucose also present).
LOG_TRANSFORM <- c("crp", "esr", "alt", "ast", "triglycerides", "glucose")

REFERENCE_LEVELS <- list(
  age_group = "Young Adult",
  sex       = "Female",
  cmv       = "Negative"
)
