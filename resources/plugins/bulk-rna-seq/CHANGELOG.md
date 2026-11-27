# Changelog — bulk-rna-seq

## 0.5.0 (2026-11-27)

- Add `rna_nine_quadrant`: cross-modal 9-quadrant integration (`Q1_rna_up_other_up`–`Q9_rna_down_other_down`) linking bulk RNA-seq DEGs with ATAC-seq/ChIP-seq/CUT&Tag/proteomics/methylation tables, with peak-to-gene aggregation (`max_abs`, `mean`, `nearest_tss`), Pearson $r$, Spearman $\rho$, 2×2 directional Fisher exact test, summary/gene tables, and labeled 3×3 scatter plot.
- Document `local/workflows/rnasplice` as the upstream alternative-splicing workflow paired with `rna_splicing`.

## 0.4.0 (2026-11-26)

- Add `rna_regulon`: `decoupler`-style Univariate Linear Model (`ULM`) and ridge-stabilized Multivariate Linear Model (`MLM`) TF/pathway regulon activity scoring from expression matrices or DE tables, with a curated 14-pathway signed footprint resource (`NFkB`, `JAK_STAT`, `TGFb`, `TNFa`, `Hypoxia`, `p53`, `MAPK`, `PI3K`, `WNT`, `Estrogen`, `Androgen`, `VEGF`, `Trail`, `EGFR`) or custom `source/target/mor` networks, analytical $t$-statistics, two-sided $p$-values, BH FDR, and activity heatmaps.
- Add `rna_splicing`: multi-format differential alternative splicing harmonization across `rMATS` (`SE`, `A5SS`, `A3SS`, `MXE`, `RI`), `SUPPA2` (`.dpsi`), and generic splicing tables with exon length extraction, reading-frame disruption classification (`frameshift_candidate` vs `in_frame_modular`), multi-hit gene aggregation, and volcano/class-breakdown plots.

## 0.3.0 (2026-11-25)

- Add `rna_coexpression`: WGCNA-style weighted gene co-expression network analysis with automatic scale-free topology $R^2$ soft-threshold ($\beta$) selection, signed/unsigned Topological Overlap Matrix (TOM) hierarchical clustering, correlated Module Eigengene merging, intramodular `kME` hub-gene ranking, and Module Eigengene–trait association testing with BH FDR.
- Add `rna_tme_deconv`: ESTIMATE-style stromal/immune enrichment scoring (`StromalScore`, `ImmuneScore`, `ESTIMATEScore`, and inferred `TumorPurity`) paired with non-negative least squares (`NNLS`) relative immune cell-type deconvolution across 10 canonical leukocyte populations (or custom reference signature matrices) and optional group differential infiltration tests.

## 0.2.0 (2026-11-25)

- Add `rna_id_map`: deterministic Ensembl/RefSeq version stripping (`ENSG00000141510.18` -> `ENSG00000141510`), GTF/GFF3 or custom 2-column gene ID mapping, optional GTF `gene_biotype` filtering, and collision-aware aggregation (`sum` for raw counts, `max_mean` for normalized matrices, `min_pvalue` for DE tables).
- Add `rna_geo_import`: parse NCBI GEO `GSE*_series_matrix.txt/.gz` files into clean `samplesheet.tsv` (parsing `!Sample_characteristics_ch1` `key: value` metadata) and `expression_matrix.tsv` with optional GPL probe-to-gene collapse and value-scale detection.
- Emit `scientificContext` metadata on `rna_id_map`, `rna_geo_import`, `rna_compare` (`valueSemantics: "de_stats"`), and `rna_pathway_score` (`valueSemantics: "zscore"`, `sampleIds`).

## 0.1.0 (2026-11-23)

Initial release: downstream bulk RNA-seq interpretation and pathway analysis plugin (`toolPrefix: rna`, agent `BulkRnaSeq`, private locked environment `plugin:rnaseq`).

- `rna_inspect`: infer matrix value kind (`raw_counts` vs TPM/FPKM/log), validate sample alignment, detect low replicates and batch-condition collinearity.
- `rna_qc`: library sizes, detected genes, `log2(CPM+1)` PCA, sample correlation heatmap, and outlier flags.
- `rna_compare`: multi-contrast DEG integration (`shared_up`, `shared_down`, `discordant`, pairwise Jaccard and `log2FC` correlation).
- `rna_ora`: offline hypergeometric ORA with explicit expressed background universe, BH FDR, and dotplot.
- `rna_gsea`: offline pre-ranked GSEA (weighted KS ES, gene-set permutation NES/p/FDR, leading-edge genes, running-sum curves).
- `rna_pathway_score`: single-sample mean gene z-score pathway activity matrix, group tests, and heatmap.
