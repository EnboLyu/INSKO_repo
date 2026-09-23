#!/usr/bin/env Rscript

# EV siINS vs EV siCONT
source("4.GSEA.R")
run_insko_gsea(
  name = "EV_siINS_vs_EV_siCONT",
  weights = c(si_siINS_vs_siCONT = 1),
  min_size = 10L,
  max_size = 500L
)
