# Changelog — proteomics

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
