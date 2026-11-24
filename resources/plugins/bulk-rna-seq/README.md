# bulk-rna-seq · 0.1.0

Downstream bulk RNA-seq interpretation and pathway analysis plugin (`toolPrefix: rna`, agent `BulkRnaSeq`) in a lightweight Python environment (`plugin:rnaseq`: Python 3.13 + `pandas/numpy/scipy/matplotlib-base`).

## Division of work with wrappers

- **Upstream alignment & quantification**: `nf-core/rna-seq` wrapper.
- **Core differential expression & time-series modeling**: `local/modules/differential-expression-{qc,deseq2,edger,limma,timeseries}` and `visualization/{pca,heatmap,volcano}` wrappers.
- **Downstream interpretation & pathway layer (this plugin, 6 `rna_*` tools)**:
  - `rna_inspect`: infer matrix value kind (`raw_counts` vs TPM/FPKM/log), validate sample alignment, detect low replicates and batch-condition collinearity.
  - `rna_qc`: library sizes, detected genes, `log2(CPM+1)` PCA, sample correlation heatmap, outlier flags.
  - `rna_compare`: multi-contrast DEG integration (`shared_up`, `shared_down`, `discordant`, pairwise Jaccard & `log2FC` correlation).
  - `rna_ora`: offline hypergeometric ORA with explicit expressed background universe, BH FDR, and dotplot.
  - `rna_gsea`: offline pre-ranked GSEA (weighted KS ES, gene-set permutation NES/p/FDR, leading-edge genes, running-sum curves).
  - `rna_pathway_score`: single-sample mean gene z-score pathway activity matrix, group tests, and heatmap.
