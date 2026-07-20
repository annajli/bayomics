# config/L3a_clinical_full_olink.R
# Bootstrap pipeline version — node maps and exclusions identical to bn_cv.
# L3a: Combined — full clinical panel + Olink proteomics
# 3 roots + 34 clinical continuous + 34 Olink continuous = 71 nodes total
#
# Built by merging the current (post-fix) L1a_clinical_full_v2.R and
# L2_olink_v2.R configs. Both source configs' exclusion decisions are
# preserved as-is here — this file does not re-derive them, it inherits them:
#   - Clinical: 8 EXCLUDE_CANDIDATES (missingness: egfr_aa/egfr_non_aa,
#     anti_ccp3/anti_ccp31, rf_iga/rf_igm results, bmi, height) + 10
#     DETERMINISTIC_EXCLUDE (mch/mchc/mcv, 5 perc_* columns, non_hdl,
#     hdl_ratio -- see L1a_clinical_full_v2.R for full justification)
#   - Olink: 1 EXCLUDE_CANDIDATE (MMP9, missingness)
# No new exclusions introduced by combining the two modalities -- verified
# no alias name collisions between the 34 clinical and 34 Olink node names
# (clinical uses lowercase snake_case, Olink uses uppercase gene symbols).
#
# NOT YET CHECKED: cross-modality collinearity (e.g. a clinical marker that's
# highly correlated with a specific Olink protein for a known biological
# reason). The within-modality checks done for L1a/L2 individually don't
# cover relationships *between* the two panels -- if this network throws
# cycle warnings, that's the most likely new source, distinct from anything
# already fixed in either source config.

NETWORK_ID <- "L3a_clinical_full_olink"
SEED <- 47L

EXCLUDE_CANDIDATES <- c(
  # -- from L1a_clinical_full_v2.R --
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

DETERMINISTIC_EXCLUDE <- c(
  # -- from L1a_clinical_full_v2.R --
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
  # -- from L2_olink_v2.R -- (none)
)

ROOT_MAP <- c(
  age_group = "subject.ageGroup",
  sex       = "subject.biologicalSex",
  cmv       = "cmv.igg_serology_interpretation"
)

# Continuous nodes (68 = 34 clinical + 34 Olink)
CONTINUOUS_MAP <- c(
  # === Clinical (34) === (from L1a_clinical_full_v2.R)
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
  triglycerides    = "clinical.lip.triglycerides",
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

# Log transform: clinical's right-skewed variables only. Olink NPX values
# are already log2-scale (no transform needed, per L2_olink_v2.R).
LOG_TRANSFORM <- c("crp", "esr", "alt", "ast", "triglycerides", "glucose")

REFERENCE_LEVELS <- list(
  age_group = "Young Adult",
  sex       = "Female",
  cmv       = "Negative"
)
