# Changelog — clinical

## 0.1.0 (2026-11-26)

Initial unified release consolidating survival biostatistics and prediction/causal/meta-analysis workflows (`toolPrefix: clin`, agent `Clinical`, private locked environment `plugin:clinical`):

- **Survival workflow (`survival-workflow`, 6 tools)**:
  - `clin_inspect`: clinical survival cohort QC, event-coding normalization, reverse-KM follow-up, EPV budget, and molecular matrix sample overlap.
  - `clin_table1`: stratified Table 1 with Shapiro-Wilk normality routing, Welch/ANOVA/Mann-Whitney/Kruskal-Wallis, Fisher/$\chi^2$, and SMD.
  - `clin_km`: Kaplan-Meier curves with log-log Greenwood 95% CIs, exact k-group log-rank test, median survival CIs, landmark rates, and Number-at-Risk subplot.
  - `clin_cox`: Breslow univariable/multivariable Cox PH regression, Harrell's C-index, Schoenfeld PH diagnostics, and HR forest plot.
  - `clin_time_roc`: leak-free K-fold L2-Cox risk scoring and IPCW time-dependent ROC evaluation.
  - `clin_competing_risk`: Aalen-Johansen CIF, naive `1 - KM` overestimation audit, Gray's test, and Fine-Gray IPCW `sHR` regression.
- **Prediction, causal & meta-analysis workflow (`prediction-causal-workflow`, 5 tools)**:
  - `clin_diagnostic_roc`: DeLong ROC variance & paired DeLong comparison, PR-AUC bootstrap CIs, and 2×2 diagnostic metrics with Wilson CIs.
  - `clin_calibration_dca`: TRIPOD calibration slope/CITL/Brier/IPA/HL diagnostics and Vickers Decision Curve Analysis (DCA).
  - `clin_psm_iptw`: positivity check, 1:$k$ caliper PSM, stabilized IPTW with ESS, Love-plot SMD audit, and causal treatment effect estimation.
  - `clin_rcs`: Harrell's Restricted Cubic Splines for linear/logistic/Cox models with overall & non-linearity Wald tests.
  - `clin_meta`: Fixed/DL/PM/HKSJ random-effects meta-analysis, 95% prediction interval, Egger's test, subgroups, Leave-One-Out, and forest/funnel plots.
- All file-producing tools emit contract-1.0.0 `.phi-artifact.json` descriptors and validate against JSON output schemas.
