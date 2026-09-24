#!/usr/bin/env Rscript

# KO siINS vs KO siCONT
source("4.0_GSEA.R")
run_insko_gsea(
  name = "KO_siINS_vs_KO_siCONT",
  weights = c(si_siINS_vs_siCONT = 1, backgroundKO.sisiINS = 1),
  min_size = 10L,
  max_size = 500L
)
