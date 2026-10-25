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
