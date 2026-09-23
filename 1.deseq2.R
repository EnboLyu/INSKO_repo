# 1. Build the DESeq2 dataset from salmon gene counts.
source("repro_helpers.R")
require_packages(c("readr", "dplyr", "DESeq2", "tibble", "apeglm"), "1.deseq2.R")

library(readr)
library(dplyr)
library(DESeq2)
library(tibble)
library(apeglm)

# ---- paths ----
counts_path <- "counts/salmon.merged.gene_counts.tsv"

out_dir <- "results/deseq2"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

dds_path     <- file.path(out_dir, "dds.rds")
cts_path     <- file.path(out_dir, "cts.rds")
coldata_path <- file.path(out_dir, "coldata.rds")
mapping_path <- "gene_name_mapping.csv"

cts_tbl <- read_tsv(counts_path)

# gene_id -> gene_name lookup, saved before we drop the column.
gene_name <- cts_tbl %>%
  select(gene_id, gene_name)
write_csv(gene_name, mapping_path)

# count matrix: gene_id as rownames, integer counts.
cts <- cts_tbl %>%
  select(-gene_name) %>%
  column_to_rownames("gene_id") %>%
  as.matrix()
mode(cts) <- "numeric"
cts <- round(cts)

coldata <- data.frame(sample = colnames(cts)) %>%
  mutate(
    background = case_when(
      grepl("KO", sample) ~ "KO",
      TRUE                ~ "EV"
    ),
    si = case_when(
      grepl("siINS", sample)  ~ "siINS",
      grepl("siCONT", sample) ~ "siCONT",
      TRUE                    ~ "Base"
    )
  )

coldata$background <- factor(coldata$background, levels = c("EV", "KO"))
coldata$si         <- factor(coldata$si, levels = c("siCONT", "siINS", "Base"))
rownames(coldata)  <- coldata$sample
coldata$sample     <- NULL


# fit model
dds <- DESeqDataSetFromMatrix(
  countData = cts,
  colData   = coldata,
  design    = ~ background * si
)
dds <- DESeq(dds)

resultsNames(dds)
# [1] "Intercept"            "background_KO_vs_EV"  "si_siINS_vs_siCONT"   "si_Base_vs_siCONT"    "backgroundKO.sisiINS"
# [6] "backgroundKO.siBase"

# QC plots
plot_dir <- file.path(out_dir, "1_deseq2")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

pdf(file.path(plot_dir, "dispersion_estimates.pdf"))
plotDispEsts(dds)
dev.off()

saveRDS(dds, dds_path)
saveRDS(cts, cts_path)
saveRDS(coldata, coldata_path)
