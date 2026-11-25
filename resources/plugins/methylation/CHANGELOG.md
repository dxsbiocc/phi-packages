# Changelog — methylation

## 0.1.0 (2026-11-25)

- Initial release of the `methylation` plugin (`toolPrefix: meth`, agent `Methylation`, skill `methylation-workflow`) in locked Python 3.13 environment `plugin:methylation` (`darwin-arm64`, `darwin-x64`, `linux-x64`), providing the downstream statistical counterpart to `nf-core/methylseq`, `bismark`, and `methyldackel` wrappers.
- Added 5 deterministic tools (`meth_inspect`, `meth_merge`, `meth_qc`, `meth_dmr`, `meth_annotate`) with JSON output schemas, `scientificContext` annotations, and contract-1.0.0 `.phi-artifact.json` descriptors.
