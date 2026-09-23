#!/usr/bin/env Rscript

# UT KO vs UT EV
source("4.GSEA.R")
run_insko_gsea(
  name = "UT_KO_vs_EV",
  weights = c(background_KO_vs_EV = 1, backgroundKO.siBase = 1),
  min_size = 15L,
  max_size = 2000L
)
