# clinical-survival · 0.1.0

Clinical cohort biostatistics and survival prognostic modeling plugin (`toolPrefix: clin`, agent `ClinicalSurvival`, private pure-conda environment `plugin:clinical`).

## Scope & division of work

- **Upstream omics quantification & DE**: `bulk-rna-seq`, `proteomics`, `genomics`, `metabolomics`, `single-cell`, or `nf-core` wrappers produce sample-level molecular scores, mutations, or expression matrices.
- **Public cancer cohort retrieval**: `cancer-omics` and `cbioportal` MCP connectors fetch TCGA / GDC / cBioPortal clinical and molecular tables.
- **This plugin (`clin_*`, 5 deterministic tools)**:
  1. `clin_inspect` — validate clinical survival tables (time positivity, inferred time unit, event coding `0/1` or `1:DECEASED`/`Dead/Alive`), reverse-KM median follow-up, event rate, Events-Per-Variable (EPV) capacity, missingness, and sample overlap with an optional molecular matrix.
  2. `clin_table1` — publication-ready stratified Baseline Characteristics Table (Table 1): Shapiro-Wilk normality check (`Mean ± SD` vs `Median [IQR]`), Welch t-test / ANOVA vs Mann-Whitney / Kruskal-Wallis, Fisher's exact vs Pearson $\chi^2$, and Standardized Mean Difference (SMD).
  3. `clin_km` — exact Kaplan-Meier product-limit survival curves with complementary log-log Greenwood 95% CIs, 2-group & multi-group exact log-rank tests, Brookmeyer-Crowley median survival CIs, landmark survival rates, and aligned Number-at-Risk subplot.
  4. `clin_cox` — Breslow partial-likelihood univariable and multivariable Cox proportional-hazards regression (`HR`, 95% CI, Wald `p`, BH FDR), Harrell's C-index, Schoenfeld residuals proportional-hazards assumption test, EPV guardrails, and publication-ready HR forest plot.
  5. `clin_time_roc` — IPCW (Inverse Probability of Censoring Weighting) time-dependent ROC curves and landmark AUCs (e.g. 1/3/5-year) for pre-computed risk scores or leak-free cross-validated L2-penalized multivariable Cox signatures.
