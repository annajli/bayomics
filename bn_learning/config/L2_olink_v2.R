# config/L2_olink.R
# L2: Olink-only — 3 roots + 35 Olink proteins
# NPX values are already log2-scale — no additional log transform needed

NETWORK_ID <- "L2_olink"
SEED <- 44L

EXCLUDE_CANDIDATES <- c("olink.MMP9")  # No olink columns are exclude_candidates
DETERMINISTIC_EXCLUDE <- character(0)

ROOT_MAP <- c(
  age_group = "subject.ageGroup",
  sex       = "subject.biologicalSex",
  cmv       = "cmv.igg_serology_interpretation"
)

CONTINUOUS_MAP <- c(
  # Inflammatory cytokines / chemokines
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
  # Chemokines (monocyte/macrophage)
  CCL2    = "olink.CCL2",
  CCL3    = "olink.CCL3",
  CCL4    = "olink.CCL4",
  CCL8    = "olink.CCL8",
  CCL11   = "olink.CCL11",
  # Aging / senescence markers
  GDF15   = "olink.GDF15",
  FGF21   = "olink.FGF21",
  CHI3L1  = "olink.CHI3L1",
  IGFBP7  = "olink.IGFBP7",
  # TNF receptor family
  TNFRSF1A = "olink.TNFRSF1A",
  TNFRSF1B = "olink.TNFRSF1B",
  # Cytotoxic / NK markers
  CD8A    = "olink.CD8A",
  GZMA    = "olink.GZMA",
  GZMB    = "olink.GZMB",
  KLRD1   = "olink.KLRD1",
  # Matrix metalloproteinases
  MMP12   = "olink.MMP12",
  # Angiogenesis
  VEGFA   = "olink.VEGFA",
  # Immune checkpoints
  PDCD1   = "olink.PDCD1",
  LAG3    = "olink.LAG3",
  HAVCR2  = "olink.HAVCR2",
  CD274   = "olink.CD274",
  # Adipokine
  LEP     = "olink.LEP"
)

NODE_MAP <- c(ROOT_MAP, CONTINUOUS_MAP)
ROOTS <- names(ROOT_MAP)

# NPX is already log2-normalized — no additional log transforms
LOG_TRANSFORM <- character(0)

REFERENCE_LEVELS <- list(
  age_group = "Young Adult",
  sex       = "Female",
  cmv       = "Negative"
)
