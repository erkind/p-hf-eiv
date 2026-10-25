# HM estimation of factor models for HF returns with EIV
#
# Minimal R script for OLS and HM regressions



# Data
R_mat <- read.csv(file.path("https://raw.githubusercontent.com/",
                            "erkind/p-hf-eiv/main/R_mat.csv"),
                  header = TRUE, sep = ",")
F_mat <- read.csv(file.path("https://raw.githubusercontent.com/",
                            "erkind/p-hf-eiv/main/F_mat.csv"),
                  header = TRUE, sep = ",")

# Model specifications
#
spec_map <- list(
  CAPM       = 1,
  FF_3F      = c(1, 3, 4),
  FFC_4F     = c(1, 3:5),
  FFC_4F_EMB = c(1, 3:5, 9),
  FFC_4F_CRD = c(1, 3:7),
  FFC_4F_OPT = c(1, 3:5, 12, 13, 16)
)

# Select specification
spec_chosen <- "FFC_4F"

# Update F_mat
F_mat <- F_mat[, spec_map[[spec_chosen]], drop = FALSE]


