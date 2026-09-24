#!/usr/bin/env Rscript

# (KO siINS response) minus (EV siINS response)
source("4.0_GSEA.R")
run_insko_gsea(
  name = "interaction_siINS_KO_minus_EV",
  weights = c(backgroundKO.sisiINS = 1),
  min_size = 10L,
  max_size = 500L
)
