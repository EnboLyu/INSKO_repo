#!/usr/bin/env Rscript

# five preranked GSEA contrasts
# Rank all genes with finite DESeq2 Wald statistics
source("repro_helpers.R")
require_packages(c("DESeq2", "fgsea", "msigdbr"), "4.0_GSEA.R")

suppressPackageStartupMessages(library(DESeq2))

gsea_dds_path <- "results/deseq2/dds.rds"
gsea_collections <- c("GO", "KEGG", "REACTOME", "HALLMARK")

gsea_coefficient_vector <- function(dds, weights) {
  coefficients <- resultsNames(dds)
  missing <- setdiff(names(weights), coefficients)
  if (length(missing) || anyDuplicated(names(weights)) ||
      !is.numeric(weights) || !length(weights)) {
    stop("Invalid DESeq2 contrast coefficients: ", paste(missing, collapse = ", "))
  }
  stats::setNames(vapply(coefficients, function(coefficient) {
    if (coefficient %in% names(weights)) weights[[coefficient]] else 0
  }, numeric(1)), coefficients)
}

gsea_rank_genes <- function(dds, weights, sample_filter = NULL) {
  contrast <- gsea_coefficient_vector(dds, weights)
  res <- results(dds, contrast = unname(contrast), alpha = 0.05)
  if (!is.null(sample_filter)) {
    selected <- which(as.character(colData(dds)$si) == sample_filter)
    if (length(selected) != 6L || !identical(sample_filter, "Base")) {
      stop("Untransfected GSEA requires the six Base samples.")
    }
    detected <- rowSums(counts(dds)[, selected, drop = FALSE]) > 0
    res <- res[detected, ]
  }
  ranked <- data.frame(gene_id = rownames(res), wald_stat = res$stat)
  ranked <- ranked[is.finite(ranked$wald_stat), , drop = FALSE]
  if (anyDuplicated(ranked$gene_id) || !nrow(ranked)) {
    stop("GSEA ranking has duplicate gene IDs or no finite statistics.")
  }
  ranked <- ranked[order(-ranked$wald_stat, ranked$gene_id), , drop = FALSE]
  rownames(ranked) <- NULL
  ranked
}

gsea_msigdb_sets <- function() {
  msig <- msigdbr::msigdbr(species = "Homo sapiens")
  gene_col <- intersect(c("ensembl_gene", "ensembl_gene_id"), names(msig))[1]
  cat_col <- intersect(c("gs_cat", "gs_collection"), names(msig))[1]
  subcat_col <- intersect(c("gs_subcat", "gs_subcollection"), names(msig))[1]
  if (anyNA(c(gene_col, cat_col, subcat_col)) || !"gs_name" %in% names(msig)) {
    stop("MSigDB columns differ from the expected msigdbr schema.")
  }
  sets <- data.frame(
    gene_id = msig[[gene_col]],
    pathway = msig$gs_name,
    category = msig[[cat_col]],
    subcategory = msig[[subcat_col]],
    stringsAsFactors = FALSE
  )
  sets <- sets[!is.na(sets$gene_id) & nzchar(sets$gene_id) &
                 !is.na(sets$pathway) & nzchar(sets$pathway), , drop = FALSE]
  sets$collection <- ifelse(sets$category == "H", "HALLMARK",
    ifelse(sets$subcategory %in% c("GO:BP", "GO:CC", "GO:MF") |
             grepl("^GOBP_|^GOCC_|^GOMF_", sets$pathway), "GO",
    ifelse(sets$subcategory == "CP:KEGG_LEGACY", "KEGG",
    ifelse(grepl("REACTOME", sets$subcategory, fixed = TRUE) |
             grepl("REACTOME", sets$pathway, fixed = TRUE), "REACTOME", NA_character_))))
  sets <- unique(sets[!is.na(sets$collection),
                      c("collection", "subcategory", "pathway", "gene_id")])
  pathways <- stats::setNames(lapply(gsea_collections, function(collection) {
    part <- sets[sets$collection == collection, , drop = FALSE]
    split(part$gene_id, part$pathway)
  }), gsea_collections)
  summary <- aggregate(gene_id ~ collection + subcategory + pathway,
                       data = sets, FUN = length)
  names(summary)[names(summary) == "gene_id"] <- "n_genes"
  summary <- summary[order(match(summary$collection, gsea_collections),
                           summary$subcategory, summary$pathway), ]
  list(pathways = pathways, summary = summary)
}

gsea_collection <- function(pathways, ranked, collection, min_size, max_size) {
  stats <- stats::setNames(ranked$wald_stat, ranked$gene_id)
  set.seed(1)
  result <- fgsea::fgsea(
    pathways = pathways, stats = stats,
    minSize = min_size, maxSize = max_size,
    nPermSimple = 10000, eps = 0
  )
  result <- as.data.frame(result)
  result$leadingEdge <- vapply(result$leadingEdge, paste, character(1), collapse = ";")
  result$collection <- collection
  result <- result[order(result$padj, -abs(result$NES), result$pathway),
                   c("collection", "pathway", "pval", "padj", "log2err",
                     "ES", "NES", "size", "leadingEdge")]
  rownames(result) <- NULL
  result
}

run_insko_gsea <- function(name, weights, min_size, max_size,
                           sample_filter = NULL, sensitivity_caps = integer()) {
  dds <- readRDS(gsea_dds_path)
  if (!inherits(dds, "DESeqDataSet") ||
      !identical(deparse(design(dds)), "~background * si")) {
    stop("GSEA requires the fitted ~ background * si DESeq2 dataset.")
  }
  ranked <- gsea_rank_genes(dds, weights, sample_filter)
  msigdb <- gsea_msigdb_sets()
  versions <- unique(msigdbr::msigdbr_collections()$db_version)
  if (length(versions) != 1L) stop("Expected one MSigDB database version.")

  out_dir <- file.path("results/gsea", name)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  rank_path <- file.path(out_dir, "ranked_genes.csv")
  summary_path <- file.path(out_dir, "msigdb_gene_set_summary.csv")
  write.csv(ranked, rank_path, row.names = FALSE)
  write.csv(msigdb$summary, summary_path, row.names = FALSE)

  results_by_collection <- lapply(gsea_collections, function(collection) {
    result <- gsea_collection(msigdb$pathways[[collection]], ranked,
                              collection, min_size, max_size)
    write.csv(result, file.path(out_dir, paste0(tolower(collection), ".csv")),
              row.names = FALSE)
    result
  })
  combined <- do.call(rbind, results_by_collection)
  rownames(combined) <- NULL
  write.csv(combined, file.path(out_dir, "all_pathways.csv"), row.names = FALSE)
  significant <- combined[!is.na(combined$padj) & combined$padj < 0.05, , drop = FALSE]
  write.csv(significant, file.path(out_dir, "significant_padj_0.05.csv"), row.names = FALSE)

  if (length(sensitivity_caps)) {
    sensitivity_by_cap <- lapply(sensitivity_caps, function(cap) {
      per_collection <- lapply(gsea_collections, function(collection) {
        result <- if (cap == max_size) {
          results_by_collection[[match(collection, gsea_collections)]]
        } else {
          gsea_collection(msigdb$pathways[[collection]], ranked,
                          collection, min_size, cap)
        }
        result[, c("collection", "pathway", "size", "NES", "padj")]
      })
      cap_result <- do.call(rbind, per_collection)
      names(cap_result)[4:5] <- paste0(c("NES_max", "padj_max"), cap)
      cap_result
    })
    sensitivity <- Reduce(function(left, right) {
      merge(left, right, by = c("collection", "pathway", "size"), all = TRUE)
    }, sensitivity_by_cap)
    write.csv(sensitivity, file.path(out_dir, "maxsize_sensitivity.csv"),
              row.names = FALSE)
  }

  metadata <- c(
    paste("contrast:", name),
    paste("coefficient_weights:", paste(names(weights), weights, sep = "=", collapse = ";")),
    paste("sample_filter:", if (is.null(sample_filter)) "none" else sample_filter),
    paste("ranked_genes:", nrow(ranked)),
    "rank_metric: DESeq2 Wald statistic, descending; gene_id breaks ties",
    paste("MSigDB_version:", versions),
    paste("msigdbr_version:", as.character(packageVersion("msigdbr"))),
    paste("fgsea_version:", as.character(packageVersion("fgsea"))),
    paste("DESeq2_version:", as.character(packageVersion("DESeq2"))),
    paste("minSize:", min_size), paste("maxSize:", max_size),
    "nPermSimple: 10000", "eps: 0", "seed_per_collection: 1",
    "FDR: BH within each contrast and collection; significant if padj < 0.05",
    paste("dds_md5:", unname(tools::md5sum(gsea_dds_path))),
    paste("rank_md5:", unname(tools::md5sum(rank_path))),
    paste("gene_set_summary_md5:", unname(tools::md5sum(summary_path))),
    "", capture.output(sessionInfo())
  )
  writeLines(metadata, file.path(out_dir, "run_metadata.txt"))
  cat(sprintf("%s: %d ranked genes, %d significant pathways\n",
              name, nrow(ranked), nrow(significant)))
  invisible(combined)
}
