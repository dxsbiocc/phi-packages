# multi-omics · 0.2.0

Cross-modality multi-omics joint integration plugin (`toolPrefix: mo`, agent
`MultiOmics`, private pure-conda environment `plugin:multiomics` with
three-platform explicit locks).

## Scope & division of work

- **Single-modality QC, normalization, and differential testing**: run first in
  the respective domain plugins (`bulk-rna-seq`, `proteomics`, `metabolomics`,
  `metagenomics`, `chromatin-profiling`, `methylation`).
- **Cross-modality integration (this plugin, 5 `mo_*` tools)**:
  1. `mo_inspect` — audit 2+ omics feature×sample matrices and optional sample
     metadata: check shared sample intersection (`n_shared_samples`), detect
     barcode mismatches (e.g., TCGA `-01A` vs patient IDs), infer per-block
     value scales (`raw_counts`, `relative_compositional`, `log_continuous`),
     and summarize shared feature IDs.
  2. `mo_rna_protein` — **Transcriptome–Proteome 9-Quadrant & Gene-Wise
     Concordance**: classify shared genes from RNA & Protein DE tables into 9
     regulatory regimes (`concordant_up/down`, `rna_only_*_buffered`,
     `protein_only_*_post_tx`, `discordant_*`), compute global log2FC
     correlation, and—when matched sample matrices are supplied—compute
     across-sample gene-wise mRNA–protein Spearman correlations with BH FDR.
  3. `mo_pair_assoc` — **Two-Block Global Coupling & Partial Correlation
     Network** (e.g., Microbiome–Metabolome or Transcriptome–Metabolome):
     modality-appropriate pre-transformation (`clr` for compositional taxa,
     `log1p`/`zscore`), **Symmetric Procrustes analysis** ($M^2$, $r$, seeded
     permutation $p$, arrow plot), **Mantel distance-matrix test**, and
     all-pairs Spearman or **covariate-adjusted Partial Spearman** network with
     BH FDR.
  4. `mo_factor_integration` — **Block-Balanced Multi-Omics Joint Latent Factor
     Decomposition** (MFA / Consensus SVD across $K \ge 2$ matched-sample
     blocks): normalizes each standardized block by its first singular value
     $\sigma_1(X_k)$ so high-dimensional blocks (20k genes) never dominate
     smaller blocks (200 metabolites), outputs shared sample factor scores,
     per-block variance explained (`shared` vs modality-specific), top signed
     loadings per block, and group association tests.
  5. `mo_subtype` — **Similarity Network Fusion (`SNF`) Multi-Omics Subtyping**:
     constructs scaled exponential sample-similarity kernels across matched
     omics blocks, iteratively fuses modality networks via K-NN message
     passing (Wang et al. 2014), selects optimal $K$ via eigengap +
     silhouette, reports per-modality NMI concordance, and renders cluster-ordered
     fused heatmaps and spectral embeddings.
