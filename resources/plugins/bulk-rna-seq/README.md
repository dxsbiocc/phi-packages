# bulk-rna-seq · 0.2.0

Downstream bulk RNA-seq interpretation, identifier/GEO bridging, and pathway analysis plugin (`toolPrefix: rna`, agent `BulkRnaSeq`) in a lightweight Python environment (`plugin:rnaseq`: Python 3.13 + `pandas/numpy/scipy/matplotlib-base`).

## Division of work with wrappers

- **Upstream alignment & quantification**: `nf-core/rna-seq` wrapper.
- **Core differential expression & time-series modeling**: `local/modules/differential-expression-{qc,deseq2,edger,limma,timeseries}` and `visualization/{pca,heatmap,volcano}` wrappers.
- **Downstream interpretation & pathway layer (this plugin, 8 `rna_*` tools)**:
  - `rna_geo_import`: parse NCBI GEO `GSE*_series_matrix.txt/.gz` files into structured `samplesheet.tsv` (`!Sample_characteristics_ch1` key-value columns) and `expression_matrix.tsv` (with optional GPL probe-to-symbol collapse).
  - `rna_id_map`: strip Ensembl/RefSeq `.version` suffixes, optionally filter by GTF `gene_biotype`, and map gene IDs via reference GTF/GFF3 or custom 2-column tables with collision-aware aggregation (`sum` / `max_mean` / `min_pvalue`).
  - `rna_inspect`: infer matrix value kind (`raw_counts` vs TPM/FPKM/log), validate sample alignment, detect low replicates and batch-condition collinearity.
  - `rna_qc`: library sizes, detected genes, `log2(CPM+1)` PCA, sample correlation heatmap, outlier flags.
  - `rna_compare`: multi-contrast DEG integration (`shared_up`, `shared_down`, `discordant`, pairwise Jaccard & `log2FC` correlation).
  - `rna_ora`: offline hypergeometric ORA with explicit expressed background universe, BH FDR, and dotplot.
  - `rna_gsea`: offline pre-ranked GSEA (weighted KS ES, gene-set permutation NES/p/FDR, leading-edge genes, running-sum curves).
  - `rna_pathway_score`: single-sample mean gene z-score pathway activity matrix, group tests, and heatmap.
