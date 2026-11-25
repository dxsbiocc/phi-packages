# Changelog — mendelian-randomization

## 0.1.0 (2026-11-25)

- Initial release of the `mendelian-randomization` plugin (`toolPrefix: mr`, agent `MendelianRandomization`, skill `mr-workflow`) in locked Python 3.13 environment `plugin:mr` (`darwin-arm64`, `darwin-x64`, `linux-x64`).
- Added 4 deterministic tools (`mr_gwas_qc`, `mr_harmonize`, `mr_estimate`, `mr_coloc`) with JSON output schemas, `scientificContext`, and `.phi-artifact.json` descriptors.
