# clinical · 0.2.0

Unified clinical biostatistics, survival modeling, prediction validation, causal inference, trial design, and meta-analysis plugin (`toolPrefix: clin`, agent `Clinical`, private pure-conda environment `plugin:clinical`).

## Scope & division of work

- **Upstream omics quantification & DE**: `bulk-rna-seq`, `proteomics`, `genomics`, `metabolomics`, `single-cell`, or `nf-core` wrappers produce sample-level molecular scores, mutations, or expression matrices.
- **Public cancer cohort retrieval**: `cancer-omics` and `cbioportal` MCP connectors fetch TCGA / GDC / cBioPortal clinical and molecular tables.
- **This plugin (`clin_*`, 14 deterministic tools across 2 skills)**:
  - **Skill 1 — `survival-workflow` (6 tools)**:
    1. `clin_inspect` — validate clinical survival tables (time positivity, inferred time unit, event coding `0/1` or `1:DECEASED`/`Dead/Alive`), reverse-KM median follow-up, event rate, Events-Per-Variable (EPV) capacity, missingness, and molecular matrix sample overlap.
    2. `clin_table1` — publication-ready stratified Baseline Characteristics Table (Table 1): Shapiro-Wilk normality check (`Mean ± SD` vs `Median [IQR]`), Welch t-test / ANOVA vs Mann-Whitney / Kruskal-Wallis, Fisher's exact vs Pearson $\chi^2$, and Standardized Mean Difference (SMD).
    3. `clin_km` — exact Kaplan-Meier product-limit survival curves with complementary log-log Greenwood 95% CIs, 2-group & multi-group exact log-rank tests, Brookmeyer-Crowley median survival CIs, landmark survival rates, and aligned Number-at-Risk subplot.
    4. `clin_cox` — Breslow partial-likelihood univariable and multivariable Cox proportional-hazards regression (`HR`, 95% CI, Wald `p`, BH FDR), Harrell's C-index, Schoenfeld residuals proportional-hazards assumption test, EPV guardrails, and publication-ready HR forest plot.
    5. `clin_time_roc` — IPCW (Inverse Probability of Censoring Weighting) time-dependent ROC curves and landmark AUCs (e.g. 1/3/5-year) for pre-computed risk scores or leak-free cross-validated L2-penalized multivariable Cox signatures.
    6. `clin_competing_risk` — Aalen-Johansen Cumulative Incidence Functions (CIF) for primary and competing events, naive `1 - KM` overestimation audit, Gray's subdistribution test, and IPCW Fine-Gray subdistribution hazard ratio (`sHR`) regression with robust sandwich CIs.
  - **Skill 2 — `prediction-causal-workflow` (8 tools)**:
    7. `clin_diagnostic_roc` — empirical ROC/PR curves, exact Sun & Xu / DeLong (1988) variance & paired DeLong $Z$-tests across biomarkers, inverse-direction detection, Youden's $J$ or fixed-cutoff 2×2 metrics with Wilson 95% CIs.
    8. `clin_calibration_dca` — TRIPOD calibration diagnostics (Brier score, IPA scaled Brier, calibration-in-the-large, calibration slope, Hosmer-Lemeshow test, $E/O$ ratio) on precomputed probabilities or leak-free K-fold OOF logistic models, plus Vickers Decision Curve Analysis (DCA) net benefit and interventions avoided.
    9. `clin_psm_iptw` — logistic propensity score positivity audit, greedy 1:$k$ caliper propensity score matching (PSM), stabilized truncated Inverse Probability of Treatment Weighting (IPTW) with Effective Sample Size (ESS), Love-plot SMD balance audit, and unadjusted/PSM/IPTW causal treatment effects (continuous, binary, survival).
    10. `clin_rcs` — Harrell's Restricted Cubic Splines (3/4/5 knots) for linear, logistic, and Cox models, overall & non-linearity Wald $\chi^2$ tests, delta-method 95% CI dose-response curves relative to `ref_value`, and inflection detection.
    11. `clin_meta` — Fixed-Effect, DerSimonian-Laird, Paule-Mandel, and HKSJ random-effects meta-analysis (precomputed effects/CIs or 2×2 counts with Haldane-Anscombe 0.5 correction), $I^2$/$Q$/$\tau^2$, 95% prediction interval, Egger's test, subgroup $Q_{\text{between}}$, Leave-One-Out sensitivity, and forest/funnel plots.
    12. `clin_sample_size` — Clinical trial & diagnostic study sample-size and power calculator across `survival_logrank` (Schoenfeld 1983), `two_proportions` (Fleiss/Casagrande superiority & Farrington-Manning 1990 non-inferiority), `two_means` (exact non-central $t$), and `diagnostic_auc` (Hanley & McNeil 1982), with dropout inflation and conservative/target/optimistic power curves.
    13. `clin_subgroup_interaction` — Subgroup forest-plot analysis with formal multiplicative interaction tests ($P_{\text{interaction}}$ via likelihood-ratio / joint Wald $\chi^2$) and Benjamini-Hochberg FDR multiplicity adjustment across Cox PH, Logistic, and Linear models.
    14. `clin_nomogram` — Multivariable Cox PH (with Breslow baseline survival at landmark horizons) or Logistic clinical prediction nomogram builder producing a 0–100 point scorecard, exact Total Points $\to$ Risk/Survival lookup table, EPV audit, and publication-ready alignment nomogram chart.
