# Step note: Cis-regulatory Peak-to-Gene linkage (`peak2gene`)

## Purpose
Link distal candidate regulatory peaks to target genes by correlating scATAC peak accessibility with scRNA expression within a configurable genomic window (`peak2gene_max_dist`, default 250 kb) of each gene TSS.

## How it works
- Pairs cells directly on shared barcodes when ATAC and RNA originate from multiome measurements, or aggregates by shared `group_by` cluster labels.
- Computes Pearson $r$ and Spearman $\rho$ for all cis peak–TSS pairs within `±max_distance`, plus a trans-chromosome empirical null distribution (`n_null`, default 200) and BH FDR correction.

## Reading the results
- `peak2gene_links.tsv`: ranked `peak`, `gene`, `tss_distance_bp`, `is_promoter`, `pearson_r`, `spearman_rho`, `pvalue`, and `fdr`.
- `figures/peak2gene-summary.png`: correlation vs TSS distance scatter and top cis links barplot.

## Pitfalls
- Correlation does not prove enhancer causality; distal enhancer–promoter links with high `pearson_r` (`fdr <= 0.05`) are candidates for CRISPRi perturbation validation.
