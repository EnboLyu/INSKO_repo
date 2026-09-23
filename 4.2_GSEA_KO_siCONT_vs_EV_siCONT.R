#!/usr/bin/env Rscript

# KO siCONT vs EV siCONT
source("4.GSEA.R")
run_insko_gsea(
  name = "KO_siCONT_vs_EV_siCONT",
  weights = c(background_KO_vs_EV = 1),
  min_size = 10L,
  max_size = 500L
)
