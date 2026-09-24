#!/usr/bin/env Rscript

input_dir <- "results/deseq2/2_shrunken"
out_dir <- "results/clustering/fc125"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

common_path <- file.path(input_dir,
  "common_significant_genes_three_first_level_comparisons.csv")
lfc_path <- file.path(input_dir, "clustering_input_shrunken_lfc.csv")
contrast_paths <- c(
  EV = file.path(input_dir, "EVsiINS_vs_EVsiCONT_all_genes_shrunken.csv"),
  KO = file.path(input_dir, "KOsiINS_vs_KOsiCONT_all_genes_shrunken.csv"),
  Base = file.path(input_dir, "KOsiCONT_vs_EVsiCONT_all_genes_shrunken.csv")
)

common <- read.csv(common_path, stringsAsFactors = FALSE)
lfc <- read.csv(lfc_path, stringsAsFactors = FALSE)
required_lfc <- c("gene_id", "gene_name", "KOsiCONT_vs_EVsiCONT",
                  "EVsiINS_vs_EVsiCONT", "KOsiINS_vs_KOsiCONT")
if (!all(c("gene_id", "gene_name") %in% names(common)) ||
    !all(required_lfc %in% names(lfc)) ||
    nrow(common) != 136L || anyNA(common$gene_id) ||
    anyDuplicated(common$gene_id) || anyDuplicated(lfc$gene_id)) {
  stop("Expected 136 unique genes and the three shrunken LFC columns.")
}

index <- match(common$gene_id, lfc$gene_id)
if (anyNA(index)) stop("A common gene is missing from the shrunken LFC matrix.")
lfc <- lfc[index, , drop = FALSE]
if (!identical(common$gene_name, lfc$gene_name)) {
  stop("Gene names differ between the common-gene and LFC tables.")
}

read_padj <- function(path, gene_ids) {
  result <- read.csv(path, stringsAsFactors = FALSE)
  if (!all(c("gene_id", "padj") %in% names(result)) ||
      anyDuplicated(result$gene_id)) {
    stop("Invalid DESeq2 results table: ", path)
  }
  values <- result$padj[match(gene_ids, result$gene_id)]
  if (anyNA(values) || any(values >= 0.05)) {
    stop("Common genes must have padj < 0.05 in every contrast: ", path)
  }
  values
}

padj <- lapply(contrast_paths, read_padj, gene_ids = common$gene_id)
cutoff <- log2(1.25)
classify <- function(values) {
  ifelse(values > cutoff, "up", ifelse(values < -cutoff, "down", "stable"))
}

genes <- data.frame(
  gene_id = common$gene_id,
  gene_name = common$gene_name,
  lfc_KOvsEV = lfc$KOsiCONT_vs_EVsiCONT,
  lfc_EV_siINS = lfc$EVsiINS_vs_EVsiCONT,
  lfc_KO_siINS = lfc$KOsiINS_vs_KOsiCONT,
  EV_padj = padj$EV,
  KO_padj = padj$KO,
  Base_padj = padj$Base,
  stringsAsFactors = FALSE
)
# if (any(!is.finite(as.matrix(genes[, c("lfc_KOvsEV", "lfc_EV_siINS",
#                                         "lfc_KO_siINS")])))) {
#   stop("The 136 genes must have finite shrunken LFCs.")
# }

genes$max_padj <- pmax(genes$EV_padj, genes$KO_padj, genes$Base_padj)
genes$dir_KOvsEV <- classify(genes$lfc_KOvsEV)
genes$dir_EV_siINS <- classify(genes$lfc_EV_siINS)
genes$dir_KO_siINS <- classify(genes$lfc_KO_siINS)
genes$direction <- paste(genes$dir_KOvsEV, genes$dir_EV_siINS,
                         genes$dir_KO_siINS, sep = "_")

groups <- split(genes, genes$direction)
summary <- do.call(rbind, lapply(names(groups), function(direction) {
  group <- groups[[direction]]
  data.frame(
    direction = direction,
    n_genes = nrow(group),
    mean_lfc_KOvsEV = mean(group$lfc_KOvsEV),
    mean_lfc_EV_siINS = mean(group$lfc_EV_siINS),
    mean_lfc_KO_siINS = mean(group$lfc_KO_siINS),
    genes = paste(sort(group$gene_name), collapse = ", "),
    stringsAsFactors = FALSE
  )
}))
summary <- summary[order(-summary$n_genes, summary$direction), , drop = FALSE]
summary$direction_cluster <- seq_len(nrow(summary))
summary <- summary[, c("direction_cluster", "direction", "n_genes",
                       "mean_lfc_KOvsEV", "mean_lfc_EV_siINS",
                       "mean_lfc_KO_siINS", "genes")]
rownames(summary) <- NULL

genes$direction_cluster <- summary$direction_cluster[
  match(genes$direction, summary$direction)]
genes <- genes[order(genes$direction_cluster, genes$max_padj), , drop = FALSE]
rownames(genes) <- NULL

write.csv(genes, file.path(out_dir, "common136_fc125_clusters.csv"),
          row.names = FALSE)
write.csv(summary, file.path(out_dir, "common136_fc125_group_summary.csv"),
          row.names = FALSE)

cat(sprintf("Grouped %d genes into %d fc125 direction groups.\n",
            nrow(genes), nrow(summary)))
