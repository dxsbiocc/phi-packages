# Changelog — clinical-survival

## 0.1.0 (2026-11-25)

Initial release: clinical biostatistics and survival prognostic modeling plugin (`toolPrefix: clin`, agent `ClinicalSurvival`, private locked environment `plugin:clinical`).

- `clin_inspect`: clinical survival table validation, event-coding normalization (`0/1`, `Dead/Alive`, `1:DECEASED`), reverse-KM median follow-up, Events-Per-Variable (EPV) capacity check, and optional molecular matrix sample overlap.
- `clin_table1`: stratified baseline characteristics Table 1 with Shapiro-Wilk normality routing (`Mean ± SD` vs `Median [IQR]`), Welch/ANOVA/Mann-Whitney/Kruskal-Wallis, Fisher's exact / $\chi^2$, and Standardized Mean Difference (SMD).
- `clin_km`: Kaplan-Meier survival curves with log-log Greenwood 95% CIs, exact k-group log-rank test, median survival 95% CIs, landmark rates, and aligned Number-at-Risk subplot.
- `clin_cox`: Breslow partial-likelihood univariable and multivariable Cox proportional-hazards regression, Harrell's C-index, Schoenfeld residual proportional-hazards assumption diagnostics, and HR forest plot.
- `clin_time_roc`: leak-free stratified K-fold L2-penalized Cox risk scoring and IPCW time-dependent ROC / landmark AUC evaluation.
- All file-producing tools emit contract-1.0.0 `.phi-artifact.json` descriptors and validate against JSON output schemas.
