# Changelog — multi-omics

## 0.1.1 (2026-11-25)

- Aligned `agents/MultiOmics.md` YAML frontmatter with the Phi Agent contract (`name`, `tools`, `skills`, `delegationMode`, `delegation`, `fallback`).

## 0.1.0 (2026-11-25)

- Initial release of the `multi-omics` plugin (`toolPrefix: mo`, agent `MultiOmics`, skill `multiomics-workflow`) in locked Python 3.13 environment `plugin:multiomics` (`darwin-arm64`, `darwin-x64`, `linux-x64`).
- Added 4 deterministic cross-omics tools with JSON output schemas, `scientificContext` annotations, and contract-1.0.0 `.phi-artifact.json` descriptors:
  - `mo_inspect`: multi-block sample intersection, barcode prefix diagnostics, scale detection, and feature overlap audit.
  - `mo_rna_protein`: 9-quadrant RNA–protein DE regulatory classification (`concordant`, `buffered`, `post_tx`, `discordant`) and matched-sample gene-wise correlation + BH FDR.
  - `mo_pair_assoc`: CLR/log1p/z-score pre-transforms, Symmetric Procrustes ($M^2$ + permutation $p$ + arrow plot), Mantel test, and covariate-adjusted Partial Spearman cross-block network.
  - `mo_factor_integration`: singular-value-balanced multi-block joint latent factor decomposition (MFA/Consensus SVD) with per-block variance explained, top loadings, and group tests.
