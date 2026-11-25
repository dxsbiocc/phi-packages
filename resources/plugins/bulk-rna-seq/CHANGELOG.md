# Changelog — bulk-rna-seq

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
