library(DESeq2)
library(dplyr)
library(tibble)
library(readr)
library(ashr)

out_dir <- "results/deseq2/2_shrunken"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# gene_id -> gene_name lookup
gene_name <- read_csv("gene_name_mapping.csv", show_col_types = FALSE) %>%
  select(gene_id, gene_name) %>%
  distinct()

dds_path <- "results/deseq2/dds.rds"
dds <- readRDS(dds_path)

required_coef <- c(
  "background_KO_vs_EV",
  "si_siINS_vs_siCONT",
  "si_Base_vs_siCONT",
  "backgroundKO.sisiINS",
  "backgroundKO.siBase"
)

missing_coef <- setdiff(required_coef, resultsNames(dds))

if (length(missing_coef) > 0) {
  stop("Missing expected DESeq2 coefficient(s): ", paste(missing_coef, collapse = ", "))
}

print(resultsNames(dds))

# Use raw DESeq2 Wald results for padj,
# then ashr-shrink the LFC for downstream effect-size / clustering use.
# Keep alpha = 0.1 for the original independent filtering; select the
# 136-gene set below using padj < 0.05 in all three contrasts.
shrink_result <- function(dds, comparison_name, coef = NULL, contrast = NULL, alpha = 0.1) {
  if (!is.null(coef) && !is.null(contrast)) {
    stop("Provide either coef or contrast, not both.")
  }
  
  if (is.null(coef) && is.null(contrast)) {
    stop("Provide either coef or contrast.")
  }
  
  if (!is.null(coef)) {
    res_raw <- results(dds, name = coef, alpha = alpha)
    res_shrunk <- lfcShrink(
      dds,
      coef = coef,
      type = "ashr",
      res = res_raw,
      quiet = TRUE
    )
  } else {
    res_raw <- results(dds, contrast = contrast, alpha = alpha)
    res_shrunk <- lfcShrink(
      dds,
      contrast = contrast,
      type = "ashr",
      res = res_raw,
      quiet = TRUE
    )
  }
  
  raw_df <- as.data.frame(res_raw) %>%
    rownames_to_column("gene_id") %>%
    transmute(
      gene_id,
      raw_log2FoldChange = log2FoldChange,
      raw_lfcSE = lfcSE,
      stat,
      pvalue,
      padj
    )
  
  shrunk_df <- as.data.frame(res_shrunk) %>%
    rownames_to_column("gene_id") %>%
    transmute(
      gene_id,
      shrunken_log2FoldChange = log2FoldChange,
      shrunken_lfcSE = lfcSE
    )
  
  out_df <- raw_df %>%
    left_join(shrunk_df, by = "gene_id") %>%
    left_join(gene_name, by = "gene_id") %>%
    relocate(gene_name, .after = gene_id) %>%
    mutate(comparison = comparison_name) %>%
    relocate(comparison, .after = gene_name)
  
  return(out_df)
}

# 1. EV siINS vs EV siCONT
res_EVsiINS_vs_EVsiCONT_df_shrink <- shrink_result(
  dds = dds,
  comparison_name = "EVsiINS_vs_EVsiCONT",
  coef = "si_siINS_vs_siCONT",
  alpha = 0.1
)

significant_EVsiINS_vs_EVsiCONT_shrink <- res_EVsiINS_vs_EVsiCONT_df_shrink %>%
  filter(!is.na(padj), padj < 0.05)

nrow(significant_EVsiINS_vs_EVsiCONT_shrink)

# 2. KO siINS vs KO siCONT
res_KOsiINS_vs_KOsiCONT_df_shrink <- shrink_result(
  dds = dds,
  comparison_name = "KOsiINS_vs_KOsiCONT",
  contrast = list(c("si_siINS_vs_siCONT", "backgroundKO.sisiINS")),
  alpha = 0.1
)

significant_KOsiINS_vs_KOsiCONT_shrink <- res_KOsiINS_vs_KOsiCONT_df_shrink %>%
  filter(!is.na(padj), padj < 0.05)

nrow(significant_KOsiINS_vs_KOsiCONT_shrink)

# 3. KO siCONT vs EV siCONT
res_KOsiCONT_vs_EVsiCONT_df_shrink <- shrink_result(
  dds = dds,
  comparison_name = "KOsiCONT_vs_EVsiCONT",
  coef = "background_KO_vs_EV",
  alpha = 0.1
)

significant_KOsiCONT_vs_EVsiCONT_shrink <- res_KOsiCONT_vs_EVsiCONT_df_shrink %>%
  filter(!is.na(padj), padj < 0.05)

nrow(significant_KOsiCONT_vs_EVsiCONT_shrink)

# Save full results
write_csv(
  res_EVsiINS_vs_EVsiCONT_df_shrink,
  file.path(out_dir, "EVsiINS_vs_EVsiCONT_all_genes_shrunken.csv")
)

write_csv(
  res_KOsiINS_vs_KOsiCONT_df_shrink,
  file.path(out_dir, "KOsiINS_vs_KOsiCONT_all_genes_shrunken.csv")
)

write_csv(
  res_KOsiCONT_vs_EVsiCONT_df_shrink,
  file.path(out_dir, "KOsiCONT_vs_EVsiCONT_all_genes_shrunken.csv")
)

# Save significant results
write_csv(
  significant_EVsiINS_vs_EVsiCONT_shrink,
  file.path(out_dir, "EVsiINS_vs_EVsiCONT_significant_shrunken.csv")
)

write_csv(
  significant_KOsiINS_vs_KOsiCONT_shrink,
  file.path(out_dir, "KOsiINS_vs_KOsiCONT_significant_shrunken.csv")
)

write_csv(
  significant_KOsiCONT_vs_EVsiCONT_shrink,
  file.path(out_dir, "KOsiCONT_vs_EVsiCONT_significant_shrunken.csv")
)

# Find genes significant in all three first-level comparisons
common_gene_ids <- Reduce(
  intersect,
  list(
    significant_EVsiINS_vs_EVsiCONT_shrink$gene_id,
    significant_KOsiINS_vs_KOsiCONT_shrink$gene_id,
    significant_KOsiCONT_vs_EVsiCONT_shrink$gene_id
  )
)

length(common_gene_ids)

common_genes_df <- tibble(gene_id = common_gene_ids) %>%
  left_join(gene_name, by = "gene_id")

write_csv(
  common_genes_df,
  file.path(out_dir, "common_significant_genes_three_first_level_comparisons.csv")
)

# Build shrunken LFC matrix for clustering
lfc_mat_df <- tibble(gene_id = common_gene_ids) %>%
  left_join(
    res_KOsiCONT_vs_EVsiCONT_df_shrink %>%
      select(gene_id, KOsiCONT_vs_EVsiCONT = shrunken_log2FoldChange),
    by = "gene_id"
  ) %>%
  left_join(
    res_EVsiINS_vs_EVsiCONT_df_shrink %>%
      select(gene_id, EVsiINS_vs_EVsiCONT = shrunken_log2FoldChange),
    by = "gene_id"
  ) %>%
  left_join(
    res_KOsiINS_vs_KOsiCONT_df_shrink %>%
      select(gene_id, KOsiINS_vs_KOsiCONT = shrunken_log2FoldChange),
    by = "gene_id"
  ) %>%
  left_join(gene_name, by = "gene_id") %>%
  relocate(gene_name, .after = gene_id)

write_csv(
  lfc_mat_df,
  file.path(out_dir, "clustering_input_shrunken_lfc.csv")
)

# Matrix version for clustering
lfc_mat <- lfc_mat_df %>%
  select(KOsiCONT_vs_EVsiCONT, EVsiINS_vs_EVsiCONT, KOsiINS_vs_KOsiCONT) %>%
  as.matrix()

rownames(lfc_mat) <- lfc_mat_df$gene_id

saveRDS(
  lfc_mat,
  file.path(out_dir, "clustering_input_shrunken_lfc_matrix.rds")
)
