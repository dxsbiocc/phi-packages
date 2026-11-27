# Multi-Omics Integration Playbook

## Stage 1 — Cohort & Modality Audit (`mo_inspect`)

- **Call**: `mo_inspect(matrices=["rna=rna_norm.tsv", "protein=prot_log2.tsv", ...], metadata=...)`
- **Read**: `n_shared_samples`, `modalities[*].inferred_scale`, `barcode_prefix_hint`.
- **Decide**:
  - If `n_shared_samples == 0` and `barcode_prefix_hint` reports 12-char TCGA patient prefix matches (`TCGA-AB-1234-01A` vs `TCGA-AB-1234`), harmonize sample headers before running downstream integration.
  - If any block has `inferred_scale == "raw_counts"`, normalize it first in its native domain plugin (or pass `log1p`).

## Stage 2 — RNA–Protein 9-Quadrant & Gene-Wise Concordance (`mo_rna_protein`)

- **Call**: `mo_rna_protein(rna_de=..., prot_de=..., rna_matrix=..., prot_matrix=..., output_dir=...)`
- **Quadrant biology**:
  - `concordant_up` / `concordant_down`: transcriptionally driven protein abundance changes.
  - `rna_only_up_buffered` / `rna_only_down_buffered`: translational buffering or stoichiometric complex degradation (common in ribosomal and OXPHOS subunits).
  - `protein_only_up_post_tx` / `protein_only_down_post_tx`: post-transcriptional regulation, ubiquitin-proteasome turnover, or secretion dynamics without mRNA change.
  - `discordant_*`: feedback loops or temporal phase lag between mRNA burst and protein accumulation.

## Stage 3 — Two-Block Global Coupling & Partial Correlation (`mo_pair_assoc`)

- **Call**: `mo_pair_assoc(matrix_x=..., matrix_y=..., transform_x="clr", covariates=[...], output_dir=...)`
- **Read**: `procrustes.m2`, `procrustes.r`, `procrustes.pvalue`, `mantel.r`, `mantel.pvalue`, `network.n_sig_edges`.
- **Red line**: When `procrustes.pvalue > 0.05` AND `mantel.pvalue > 0.05`, state clearly that global multivariate geometry does not couple significantly across the two blocks, and treat isolated FDR-passing feature pairs as localized hypotheses requiring validation.

## Stage 4 — Joint Latent Factor Decomposition (`mo_factor_integration`)

- **Call**: `mo_factor_integration(matrices=[...], metadata=..., group_col=..., n_factors=5, output_dir=...)`
- **Principle**: Multiple Factor Analysis (MFA / Escofier-Pagès scaling):
  1. Each modality $X_k \in \mathbb{R}^{n \times p_k}$ is pre-transformed, median-imputed for residual NAs, and feature-standardized (zero mean, unit variance).
  2. Each standardized block is divided by its first singular value $\sigma_1(X_k)$ so its maximum axial inertia equals $1.0$.
  3. SVD on concatenated $[ \tilde{X}_1, \dots, \tilde{X}_K ]$ yields shared sample scores $F$ and per-block variance explained $R^2_{k,f} = \| \hat{X}_{k,f} \|_F^2 / \| \tilde{X}_k \|_F^2$, distinguishing **shared multi-omics factors** (high $R^2$ in $\ge 2$ blocks) from **modality-specific factors** (high $R^2$ in 1 block only).

## Stage 5 — Similarity Network Fusion Subtyping (`mo_subtype`)

- **Call**: `mo_subtype(matrices=["rna=...", "protein=...", ...], metadata=..., group_col=..., output_dir=...)`
- **Read**: `k_clusters`, `silhouette_score`, `modalities[*].modality_concordance_nmi`, `group_association`.
- **Decide & Red Lines**:
  - Prefer `silhouette_score >= 0.25` and check `modality_concordance_nmi`: if one modality has near-zero NMI with the fused subtypes while others align strongly, report that fused clustering is driven primarily by the concordant modalities.
