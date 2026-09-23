#!/usr/bin/env Rscript

# ORA of 136 genes
source("repro_helpers.R")
require_packages(c("DESeq2", "gprofiler2"), "3.ORA.R")

suppressPackageStartupMessages({
  library(DESeq2)
  library(gprofiler2)
})

dds_path <- "results/deseq2/dds.rds"
genes_path <- "results/deseq2/2_shrunken/common_significant_genes_three_first_level_comparisons.csv"
out_dir <- "results/ora"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

dds <- readRDS(dds_path)
contrasts <- list(
  EVsiINS_vs_EVsiCONT = results(dds, name = "si_siINS_vs_siCONT", alpha = 0.1),
  KOsiINS_vs_KOsiCONT = results(
    dds,
    contrast = list(c("si_siINS_vs_siCONT", "backgroundKO.sisiINS")),
    alpha = 0.1
  ),
  KOsiCONT_vs_EVsiCONT = results(dds, name = "background_KO_vs_EV", alpha = 0.1)
)

# The custom background includes only genes with a defined padj in every test.
testable_ids <- lapply(contrasts, function(res) {
  rownames(res)[!is.na(res$padj)]
})
background_ids <- sort(Reduce(intersect, testable_ids))

significant_ids <- lapply(contrasts, function(res) {
  rownames(res)[!is.na(res$padj) & res$padj < 0.05]
})
expected_ids <- sort(Reduce(intersect, significant_ids))

genes <- read.csv(genes_path, stringsAsFactors = FALSE)
if (!all(c("gene_id", "gene_name") %in% names(genes)) ||
    anyNA(genes$gene_id) || anyDuplicated(genes$gene_id) ||
    !setequal(genes$gene_id, expected_ids) || length(expected_ids) != 136L) {
  stop("The 136-gene CSV does not match padj < 0.05 in all three DESeq2 contrasts.")
}

genes <- genes[order(genes$gene_id), c("gene_id", "gene_name")]
write.csv(genes, file.path(out_dir, "whole136_gene_list.csv"), row.names = FALSE)
write.csv(data.frame(gene_id = background_ids),
          file.path(out_dir, "testable_background_gene_ids.csv"), row.names = FALSE)
dump_session_info(out_dir, "ora")

cat(sprintf("ORA query: %d genes; testable background: %d genes\n",
            nrow(genes), length(background_ids)))

g <- gost(
  query = genes$gene_id,
  organism = "hsapiens",
  custom_bg = background_ids,
  domain_scope = "custom",
  sources = c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC"),
  correction_method = "g_SCS",
  evcodes = TRUE,
  significant = TRUE
)

if (is.null(g) || is.null(g$result) || nrow(g$result) == 0L ||
    is.null(g$meta$version)) {
  stop("g:Profiler returned no significant results or no database version.")
}

enrichment <- g$result[order(g$result$p_value),
  c("source", "term_id", "term_name", "p_value", "intersection_size",
    "term_size", "intersection")]
write.csv(enrichment, file.path(out_dir, "enrichment_whole136.csv"), row.names = FALSE)

provenance <- c(
  paste("gprofiler_db_version:", g$meta$version),
  paste("gprofiler_query_time:", g$meta$timestamp),
  paste("gprofiler2_version:", as.character(packageVersion("gprofiler2"))),
  paste("DESeq2_version:", as.character(packageVersion("DESeq2"))),
  paste("R_version:", R.version.string),
  "DESeq2_results_alpha: 0.1",
  "gene_padj_cutoff: 0.05 in all three contrasts",
  "organism: hsapiens",
  "sources: GO:BP,GO:MF,GO:CC,KEGG,REAC",
  "correction_method: g_SCS",
  "domain_scope: custom",
  paste("submitted_query_n:", nrow(genes)),
  paste("submitted_background_n:", length(background_ids)),
  paste("effective_query_n:", paste(unique(g$result$query_size), collapse = ",")),
  paste("effective_background_n:",
        paste(unique(g$result$effective_domain_size), collapse = ",")),
  paste("dds_md5:", unname(tools::md5sum(dds_path))),
  paste("gene_list_md5:", unname(tools::md5sum(genes_path)))
)
writeLines(provenance, file.path(out_dir, "gprofiler_provenance.txt"))

cat(sprintf("Wrote %d significant ORA terms to %s\n", nrow(enrichment), out_dir))
