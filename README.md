# INSKO RNA-seq analysis

R scripts and output files for RNA-seq analysis of 18 samples of two cell backgrounds, three treatment conditions, and three biological replicates per condition.

## Input and methods

`counts/salmon.merged.gene_counts.tsv`: the gene count input.
`gene_name_mapping.csv`: gene identifiers mapping to gene names.

DESeq2 fits all samples with the model `~ background * si`, significance threshold is defined as `padj < 0.05`.
ashr is used to shrink log2foldchange estimates.

Over-representation analysis (ORA) uses g:Profiler.

Gene set enrichment analysis (GSEA) uses fgsea and gene sets from the Molecular Signatures Database (MSigDB).
GSEA ranks genes by the DESeq2 Wald test statistic.

## Scripts

Use R 4.6.0 and the package versions in `renv.lock`.
Run scripts from the repository root.

1. Run `1.deseq2.R` to fit the model.
2. Run `2.deseq2_re_ashr.R` to export comparisons and the common gene set.
3. Use `3.ORA.R` for ORA, `4.1_*.R` to `4.5_*.R` for GSEA.
4. Use `5.cluster_fc125.R` for direction groups.

## Rresults

| Folder | Contents |
| --- | --- |
| [results/deseq2/](results/deseq2/dds.rds) | Model|
| [results/deseq2/2_shrunken/](results/deseq2/2_shrunken/) | Full and significant gene tables, and clustering input. |
| [results/ora/](results/ora/) | Enrichment results, background genes. |
| [results/gsea/](results/gsea/) | GSEA results, one folder per comparison. |
| [results/clustering/fc125/](results/clustering/fc125/) | Gene grouping results. |
