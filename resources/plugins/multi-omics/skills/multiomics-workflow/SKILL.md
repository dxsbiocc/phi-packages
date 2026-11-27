---
name: multiomics-workflow
description: >-
  Cross-modality multi-omics integration workflow in Python: inspect shared
  samples, barcode formatting, and value scales across 2+ omics matrices; run
  RNA-protein 9-quadrant DE concordance (translational buffering vs
  post-transcriptional regulation) and matched-sample gene-wise correlation;
  test global two-block coupling via Symmetric Procrustes and Mantel permutation
  tests plus covariate-adjusted Partial Spearman cross-omics networks (e.g.
  microbiome-metabolome with CLR); and fit block-balanced joint latent factor
  models (MFA/Consensus SVD) reporting per-block variance explained and top
  multi-omics loadings, and Similarity Network Fusion (SNF) patient/sample
  subtyping with silhouette QC. Use for 多组学联合, multi-omics, 转录蛋白联合, 九宫格,
  9-quadrant, 菌群代谢联合, Procrustes, Mantel, 偏相关网络, MOFA, CCA, SNF, or cross-omics
  integration.
phi:
  environment: plugin:multiomics
  attachTo:
    - MultiOmics
  scripts:
    - name: inspect
      description: >-
        Inspect two or more omics feature x sample matrices (optionally labeled
        as modality=path) and optional sample metadata: report per-modality
        dimensions, missingness, value scale, exact sample intersection, barcode
        prefix match hints, and shared feature IDs. Call first before sample-level
        integration.
      run: [python, ./scripts/mo_inspect.py]
      args:
        type: object
        additionalProperties: false
        required: [matrices]
        properties:
          matrices:
            type: array
            items:
              type: string
            description: Two or more feature x sample TSV/CSV paths, optionally labeled 'modality=path'.
          metadata:
            type: string
            format: input-path
            description: Optional sample metadata TSV/CSV.
          sample_col:
            type: string
            description: Sample ID column in metadata (default first column).
          group_col:
            type: string
            description: Optional grouping column in metadata to check per-group matched sample sizes.
      approval: read
      output: ./schemas/inspect-result.json
      timeoutSeconds: 180
    - name: rna_protein
      description: >-
        Integrate transcriptome and proteome differential expression tables
        (9-quadrant regulatory classification: concordant, translational
        buffering, post-transcriptional/turnover, discordant) and/or matched
        sample expression & abundance matrices (across-sample per-gene Spearman
        correlation + BH FDR).
      run: [python, ./scripts/mo_rna_protein.py]
      args:
        type: object
        additionalProperties: false
        required: [output_dir]
        properties:
          rna_de:
            type: string
            format: input-path
            description: Bulk RNA-seq DE table TSV/CSV (gene, log2FC, padj/FDR).
          prot_de:
            type: string
            format: input-path
            description: Proteomics DE table TSV/CSV (gene/protein, log2FC, padj/FDR).
          rna_matrix:
            type: string
            format: input-path
            description: Optional gene x sample normalized RNA matrix TSV/CSV for across-sample correlation.
          prot_matrix:
            type: string
            format: input-path
            description: Optional gene/protein x sample log2 protein matrix TSV/CSV for across-sample correlation.
          gene_map:
            type: string
            format: input-path
            description: Optional 2-column ID mapping TSV/CSV mapping protein IDs to RNA gene symbols.
          rna_fc_cutoff:
            type: number
            description: Absolute RNA log2FC threshold for 9-quadrant classification (default 1.0).
          prot_fc_cutoff:
            type: number
            description: Absolute Protein log2FC threshold for 9-quadrant classification (default 0.585, i.e. 1.5-fold).
          fdr_cutoff:
            type: number
            description: Adjusted p-value / FDR threshold for both modalities (default 0.05).
          output_dir:
            type: string
            format: project-path
            description: Output directory for 9-quadrant table, gene-wise correlations, and scatter plots.
      approval: write
      output: ./schemas/rna-protein-result.json
      timeoutSeconds: 300
    - name: pair_assoc
      description: >-
        Test global coupling and feature-pair associations between two matched
        omics matrices (e.g. microbiome x metabolome or RNA x metabolome): runs
        Symmetric Procrustes analysis (M^2, correlation r, permutation p, arrow
        plot), Mantel distance-matrix permutation test, and all-pairs Spearman
        or covariate-adjusted Partial Spearman network with BH FDR.
      run: [python, ./scripts/mo_pair_assoc.py]
      args:
        type: object
        additionalProperties: false
        required: [matrix_x, matrix_y, output_dir]
        properties:
          matrix_x:
            type: string
            format: input-path
            description: First omics feature x sample matrix TSV/CSV (e.g. microbiome or RNA).
          matrix_y:
            type: string
            format: input-path
            description: Second omics feature x sample matrix TSV/CSV (e.g. metabolome or proteome).
          name_x:
            type: string
            description: Display label for matrix_x (default block_X).
          name_y:
            type: string
            description: Display label for matrix_y (default block_Y).
          transform_x:
            type: string
            enum: [auto, none, log1p, clr, zscore]
            description: Pre-transformation for matrix_x (use clr for compositional relative abundances; default auto).
          transform_y:
            type: string
            enum: [auto, none, log1p, clr, zscore]
            description: Pre-transformation for matrix_y (default auto).
          metadata:
            type: string
            format: input-path
            description: Optional sample metadata TSV/CSV for covariate adjustment or Procrustes group coloring.
          sample_col:
            type: string
            description: Sample ID column in metadata (default first column).
          group_col:
            type: string
            description: Optional group column for coloring Procrustes plots.
          covariates:
            type: array
            items:
              type: string
            description: Optional metadata columns to partial out via OLS before computing Partial Spearman correlations.
          top_x:
            type: integer
            minimum: 5
            description: Top variable features retained from matrix_x for pairwise network screening (default 100).
          top_y:
            type: integer
            minimum: 5
            description: Top variable features retained from matrix_y for pairwise network screening (default 100).
          fdr_cutoff:
            type: number
            description: BH FDR cutoff for significant cross-block edges (default 0.05).
          permutations:
            type: integer
            minimum: 99
            description: Permutations for Procrustes and Mantel p-values (default 999).
          seed:
            type: integer
            description: Random seed (default 42).
          output_dir:
            type: string
            format: project-path
            description: Output directory for Procrustes/Mantel results, edge table, and figures.
      approval: write
      output: ./schemas/pair-assoc-result.json
      timeoutSeconds: 600
    - name: factor_integration
      description: >-
        Fit a block-balanced multi-omics joint latent factor decomposition
        (MFA / Consensus SVD) across 2 or more matched-sample omics matrices:
        weights each standardized block by 1/sigma_1(X_k) so high-dimensional
        modalities do not dominate, outputs shared sample factor scores,
        per-block R^2 variance explained, top signed feature loadings per block,
        optional group association tests, and multi-panel diagnostic figures.
      run: [python, ./scripts/mo_factor_integration.py]
      args:
        type: object
        additionalProperties: false
        required: [matrices, output_dir]
        properties:
          matrices:
            type: array
            items:
              type: string
            description: Two or more feature x sample TSV/CSV paths, optionally formatted 'modality=path'.
          transforms:
            type: array
            items:
              type: string
            description: Optional per-block transforms ('auto', 'none', 'log1p', 'clr') matching matrices order.
          metadata:
            type: string
            format: input-path
            description: Optional sample metadata TSV/CSV for factor-phenotype testing and plot coloring.
          sample_col:
            type: string
            description: Sample ID column in metadata (default first column).
          group_col:
            type: string
            description: Optional grouping column to test against each latent factor (ANOVA/Kruskal + BH FDR).
          n_factors:
            type: integer
            minimum: 2
            maximum: 15
            description: Number of joint latent factors to extract (default 5).
          max_features_per_block:
            type: integer
            minimum: 10
            description: Retain up to N highest-variance features per block before SVD (default 2000).
          top_loadings:
            type: integer
            minimum: 5
            description: Top absolute-loading features exported per block and factor (default 25).
          output_dir:
            type: string
            format: project-path
            description: Output directory for factor scores, variance explained, loadings, and plots.
      approval: write
      output: ./schemas/factor-integration-result.json
      timeoutSeconds: 600
    - name: subtype
      description: >-
        Discover multi-omics molecular subtypes via Wang et al. (2014) Similarity
        Network Fusion (SNF): constructs scaled exponential affinity networks across
        2+ matched modalities, performs iterative K-NN message-passing network fusion,
        auto-selects optimal K via normalized Laplacian eigengap + silhouette (or
        uses fixed k_clusters), reports per-modality NMI concordance, and exports
        cluster assignments, fused similarity matrix, and heatmap + spectral plots.
      run: [python, ./scripts/mo_subtype.py]
      args:
        type: object
        additionalProperties: false
        required: [matrices, output_dir]
        properties:
          matrices:
            type: array
            items:
              type: string
            description: Two or more feature x sample TSV/CSV paths, optionally formatted 'modality=path'.
          output_dir:
            type: string
            format: project-path
            description: Output directory for snf_cluster_assignments.tsv, snf_fused_similarity.tsv, and snf_fused_heatmap_embedding.png.
          transforms:
            type: array
            items:
              type: string
            description: Optional per-block transforms ('auto', 'none', 'log1p', 'clr') matching matrices order.
          metadata:
            type: string
            format: input-path
            description: Optional sample metadata TSV/CSV for external label concordance checks.
          sample_col:
            type: string
            description: Sample ID column in metadata (default first column).
          group_col:
            type: string
            description: Optional reference label column in metadata to evaluate ARI, NMI, and chi-square association.
          k_clusters:
            type: integer
            minimum: 2
            maximum: 8
            description: Optional fixed number of molecular subtypes (if omitted, selected via eigengap + silhouette).
          n_neighbors:
            type: integer
            minimum: 3
            maximum: 50
            description: Local K-nearest neighbors bandwidth for SNF affinity & sparse kernel construction (default 15).
          n_iters:
            type: integer
            minimum: 1
            maximum: 50
            description: Iterative SNF cross-modality network fusion steps (default 15).
          seed:
            type: integer
            description: Deterministic random seed for spectral K-means clustering (default 42).
      approval: write
      output: ./schemas/subtype-result.json
      timeoutSeconds: 600
---

# Multi-Omics Integration Workflow

Read [references/analysis-playbook.md](references/analysis-playbook.md) for
cross-omics normalization rules, 9-quadrant interpretation, compositional CLR
handling, confounder adjustment, and block-weight balancing.

## Workflow Chain

```text
Single-omics matrices / DE tables (from bulk-rna-seq, proteomics, metabolomics, metagenomics)
  ├── mo_inspect ────────────▶ verify shared sample intersection, barcode format & value scales
  ├── mo_rna_protein ────────▶ 9-quadrant RNA-Protein DE concordance + across-sample r & FDR
  ├── mo_pair_assoc ─────────▶ Symmetric Procrustes + Mantel + Partial Spearman cross-block network
  ├── mo_factor_integration ─▶ MFA / Consensus SVD joint latent factors + per-block R² & loadings
  └── mo_subtype ────────────▶ Similarity Network Fusion (SNF) multi-omics subtypes + concordance NMI
```

## Hard Red Lines

1. **Never correlate unaligned samples**: require $\ge 6$ shared biological samples for `mo_pair_assoc` / `mo_factor_integration` ($\ge 15$ recommended for stable Partial Spearman networks).
2. **Apply `clr` to compositional microbiome profiles**: relative abundances sum to 1 (or 100), creating negative closure bias; always pass `transform_x=clr` for taxonomic relative abundances.
3. **Adjust for disease/group confounders in cross-block networks**: pass `--covariates` in `mo_pair_assoc` when case-control separation drives both modalities simultaneously.
4. **Equalize block inertia**: `mo_factor_integration` scales each feature-standardized block by $1/\sigma_1(X_k)$ so a 20,000-gene transcriptome and a 150-metabolite panel have equal opportunity to shape the first joint factor.
