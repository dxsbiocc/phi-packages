# Changelog — proteomics

## 0.2.1 (2026-11-24)

- Add contract-1.0.0 `.phi-artifact.json` descriptors and `artifacts` lists across all file-producing tools (`prot_rollup`, `prot_qc_impute`, `prot_de`, `prot_ptm`, `prot_enrichment`, `prot_network`).
- Add `scientificContext` annotations (`log_transformed` + `sampleIds` on `prot_rollup` and `prot_qc_impute`; `de_stats` on `prot_de` and `prot_ptm`) and update output JSON schemas.

## 0.2.0 (2026-11-24)

- `prot_rollup`: added peptide/precursor-to-protein rollup (`rrollup` reference-scaled log2 rollup, `top3`/Hi3, `median`, `sum`) with `unique_only` vs Occam's `razor` shared-peptide handling, contaminant removal, per-protein peptide counts, and histogram output.

## 0.1.0 (2026-11-24)

Initial release: pure conda-forge downstream quantitative MS proteomics and
PTM/phosphoproteomics plugin (`toolPrefix: prot`).

- `prot_inspect`: format auto-detection (MaxQuant, DIA-NN, FragPipe, generic),
  contaminant/reverse detection, scale inference, missingness profiling, and
  batch-condition design checks.
- `prot_qc_impute`: contaminant removal, log2 transform, group-aware valid-value
  filtering, median/quantile normalization, Perseus-style MinProb (MNAR) and
  kNN (MAR) imputation, plus missingness/boxplot/PCA QC plots.
- `prot_de`: empirical-Bayes moderated t-test (Smyth-style variance shrinkage),
  Welch t-test, and Mann-Whitney U with BH FDR, volcano and MA plots.
- `prot_ptm`: PTM/phosphosite differential analysis with optional total-protein
  normalization (`site_log2 - protein_log2`) and Kinase-Substrate Enrichment
  Analysis (KSEA).
- `prot_enrichment`: offline hypergeometric ORA with explicit detected-proteome
  background universe and deterministic pre-ranked GSEA against GMT files.
- `prot_network`: protein co-abundance correlation network, module extraction,
  hub scoring, and clustered heatmap.
