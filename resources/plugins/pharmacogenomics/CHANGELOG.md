# Changelog — pharmacogenomics

## 0.2.0 (2026-11-26)

- Add `pgx_hts_qc`: high-throughput screening (96/384-well) plate QC computing Zhang's $Z'$-factor, `SSMD`, signal-to-background (`S/B`), control `CV%`, row/column ANOVA and edge-effect spatial bias diagnostics, Tukey two-way median polish `B-score`, `% inhibition`, robust MAD `z_robust` hit calling, and plate spatial heatmaps.
- Add `pgx_predict_cohort`: cross-cohort cell-line-to-patient pharmacogenomic translation (`oncoPredict`/Geeleher-inspired) with cross-platform gene-wise Z-score/empirical-Bayes harmonization, 5-fold cross-validated `RidgeCV` modeling, out-of-fold Pearson `r` / Spearman $\rho$ / RMSE evaluation, `low_predictability_warning` flagging, and patient sensitivity overview plots.

## 0.1.0 (2026-11-25)

Initial release: pharmacogenomics, dose-response curve modeling, 2D drug synergy, and drug-sensitivity biomarker discovery plugin (`toolPrefix: pgx`, agent `Pharmacogenomics`, private locked environment `plugin:pgx`).

- `pgx_dose_response`: 4-parameter log-logistic Hill curve fitting (`EC50`, absolute `IC50`, `E0`, `Einf`, Hill slope, `AUC`, `AAC`, $R^2$) with explicit `beyond_tested_range` censoring guards when `IC50` exceeds maximum tested concentration.
- `pgx_synergy`: 2D drug-combination viability matrix modeling across **ZIP**, **Bliss**, **Loewe**, and **HSA** reference models, peak $3 \times 3$ synergy window detection, and dose-response / ZIP synergy landscape heatmaps.
- `pgx_biomarker`: genomic alteration (MAF or binary mutation/CNV matrix) and transcriptomic/pathway sensitivity biomarker screening with explicit drug-metric polarity (`higher_is_sensitive`) and Benjamini–Hochberg FDR correction.
