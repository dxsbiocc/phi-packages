---
name: survival-workflow
description: >-
  Clinical cohort biostatistics and survival prognostic workflow in Python:
  audit follow-up time, event/censoring encoding, reverse-KM follow-up, and EPV
  capacity; generate publication-ready stratified Baseline Characteristics Table
  1 (Shapiro-Wilk normality routing, Welch/ANOVA/Mann-Whitney/Kruskal-Wallis,
  Fisher/chi-square, SMD); plot Kaplan-Meier curves with log-log Greenwood 95%
  CIs, exact log-rank tests, and aligned Number-at-Risk tables; fit univariable
  and multivariable Breslow Cox proportional-hazards models with Harrell's
  C-index, Schoenfeld PH assumption tests, and HR forest plots; and evaluate
  prognostic risk scores or leak-free cross-validated L2-Cox signatures via IPCW
  time-dependent ROC (1/3/5-year AUC). Use for 生存分析, Kaplan-Meier, Log-rank,
  Cox回归, 风险比, HR森林图, 预后模型, 时间依赖ROC, Table 1, 基线特征表, clinical
  biostatistics, or survival analysis.
phi:
  environment: plugin:clinical
  attachTo:
    - ClinicalSurvival
  scripts:
    - name: inspect
      description: >-
        Inspect a clinical cohort table (TSV/CSV) and optional molecular matrix:
        normalize event coding (0/1, Dead/Alive, 1:DECEASED), check survival
        time positivity and infer units (days/months/years), estimate reverse-KM
        median follow-up, compute Events-Per-Variable (EPV) budget, and summarize
        missingness and molecular sample overlap. Call first on every cohort.
      run: [python, ./scripts/clin_inspect.py]
      args:
        type: object
        additionalProperties: false
        required: [clinical, time_col, event_col]
        properties:
          clinical:
            type: string
            format: input-path
            description: Clinical metadata table TSV/CSV.
          time_col:
            type: string
            description: Follow-up time column name (must be > 0 for valid records).
          event_col:
            type: string
            description: Event indicator column (1/Dead/DECEASED = event, 0/Alive/LIVING = censored).
          sample_col:
            type: string
            description: Sample/patient ID column (default first column).
          matrix:
            type: string
            format: input-path
            description: Optional feature x sample molecular matrix TSV/CSV to check sample overlap.
      approval: read
      output: ./schemas/inspect-result.json
      timeoutSeconds: 180
    - name: table1
      description: >-
        Build a publication-ready Baseline Characteristics Table (Table 1),
        optionally stratified by group_col: auto-detects continuous vs
        categorical variables, tests normality via Shapiro-Wilk to report both
        Mean +/- SD and Median [IQR], runs Welch/ANOVA or Mann-Whitney/Kruskal
        and Fisher's exact or chi-square tests, and computes two-group SMD.
      run: [python, ./scripts/clin_table1.py]
      args:
        type: object
        additionalProperties: false
        required: [clinical, variables, output_dir]
        properties:
          clinical:
            type: string
            format: input-path
            description: Clinical metadata table TSV/CSV.
          variables:
            type: array
            items:
              type: string
            description: Clinical or molecular columns to summarize in Table 1.
          group_col:
            type: string
            description: Optional stratification column (e.g. treatment arm, stage, risk group).
          categorical_cols:
            type: array
            items:
              type: string
            description: Optional explicit list of variables to treat as categorical.
          continuous_cols:
            type: array
            items:
              type: string
            description: Optional explicit list of variables to treat as continuous.
          output_dir:
            type: string
            format: project-path
            description: Output directory for table1.tsv and overview figure.
      approval: write
      output: ./schemas/table1-result.json
      timeoutSeconds: 300
    - name: km
      description: >-
        Compute exact Kaplan-Meier survival curves with log-log Greenwood 95%
        confidence intervals, Brookmeyer-Crowley median survival CIs, landmark
        survival probabilities, exact 2-group or multi-group log-rank test, and
        a two-panel KM curve + aligned Number-at-Risk figure. Supports either a
        categorical group_col or pre-specified quantile cutoffs of biomarker_col.
      run: [python, ./scripts/clin_km.py]
      args:
        type: object
        additionalProperties: false
        required: [clinical, time_col, event_col, output_dir]
        properties:
          clinical:
            type: string
            format: input-path
            description: Clinical metadata table TSV/CSV.
          time_col:
            type: string
            description: Follow-up time column.
          event_col:
            type: string
            description: Event indicator column.
          sample_col:
            type: string
            description: Optional sample ID column (required when joining --matrix).
          group_col:
            type: string
            description: Categorical grouping column in clinical table.
          biomarker_col:
            type: string
            description: Continuous variable (in clinical table or row in --matrix) to stratify by quantile/cutoff.
          matrix:
            type: string
            format: input-path
            description: Optional feature x sample matrix if biomarker_col is a molecular feature.
          cutoff:
            type: string
            description: Cutoff rule for biomarker_col — median (default), tertile, quartile, or an explicit numeric value.
          landmarks:
            type: array
            items:
              type: number
            description: Optional landmark follow-up times (e.g. [12, 36, 60]) for survival rate reporting.
          time_label:
            type: string
            description: Axis label for time unit (e.g. Months, Years, Days).
          output_dir:
            type: string
            format: project-path
            description: Output directory for km_summary.tsv, km_curves.tsv, and km-curve.png.
      approval: write
      output: ./schemas/km-result.json
      timeoutSeconds: 300
    - name: cox
      description: >-
        Fit Breslow partial-likelihood univariable and/or multivariable Cox
        proportional-hazards regression models: outputs Hazard Ratios (HR), 95%
        Wald CIs, p-values, BH FDR (for univariable screens), Harrell's C-index,
        Schoenfeld residual proportional-hazards diagnostics, EPV checks, and a
        publication-ready HR forest plot.
      run: [python, ./scripts/clin_cox.py]
      args:
        type: object
        additionalProperties: false
        required: [clinical, time_col, event_col, covariates, output_dir]
        properties:
          clinical:
            type: string
            format: input-path
            description: Clinical metadata table TSV/CSV.
          time_col:
            type: string
            description: Follow-up time column.
          event_col:
            type: string
            description: Event indicator column.
          sample_col:
            type: string
            description: Optional sample ID column (used when joining --matrix features).
          covariates:
            type: array
            items:
              type: string
            description: Clinical columns and/or molecular feature IDs (from --matrix) to model.
          matrix:
            type: string
            format: input-path
            description: Optional feature x sample matrix providing molecular covariates.
          mode:
            type: string
            enum: [both, multivariable, univariable]
            description: Fit univariable screen, multivariable model, or both (default both).
          standardize_continuous:
            type: boolean
            description: Z-score standardize continuous predictors so HR reflects per-1-SD change (default false).
          output_dir:
            type: string
            format: project-path
            description: Output directory for Cox tables and cox-forest.png.
      approval: write
      output: ./schemas/cox-result.json
      timeoutSeconds: 300
    - name: time_roc
      description: >-
        Evaluate a pre-computed continuous risk score (score_col) OR fit a
        leak-free stratified K-fold cross-validated L2-penalized multivariable
        Cox risk signature (feature_cols), computing IPCW time-dependent ROC
        curves and AUCs at clinical landmarks (e.g. 1, 3, 5 years or 12, 36, 60
        months) plus out-of-fold median high/low risk KM log-rank separation.
      run: [python, ./scripts/clin_time_roc.py]
      args:
        type: object
        additionalProperties: false
        required: [clinical, time_col, event_col, landmarks, output_dir]
        properties:
          clinical:
            type: string
            format: input-path
            description: Clinical metadata table TSV/CSV.
          time_col:
            type: string
            description: Follow-up time column.
          event_col:
            type: string
            description: Event indicator column.
          sample_col:
            type: string
            description: Sample ID column (default first column).
          score_col:
            type: string
            description: Existing continuous risk score column in clinical table (mutually optional with feature_cols).
          feature_cols:
            type: array
            items:
              type: string
            description: Candidate clinical/molecular features for K-fold L2-penalized Cox out-of-fold scoring.
          matrix:
            type: string
            format: input-path
            description: Optional feature x sample molecular matrix when feature_cols include genes/proteins.
          landmarks:
            type: array
            items:
              type: number
            description: Landmark follow-up times (e.g. [1, 3, 5] or [12, 36, 60]) for IPCW time-dependent ROC.
          n_folds:
            type: integer
            minimum: 3
            maximum: 10
            description: Cross-validation folds when fitting feature_cols (default 5).
          l2_penalty:
            type: number
            description: Ridge L2 penalty lambda for Cox signature fitting (default 1.0).
          seed:
            type: integer
            description: Random seed for fold assignment (default 42).
          output_dir:
            type: string
            format: project-path
            description: Output directory for landmark_auc.tsv, risk_scores.tsv, and time-dependent-roc.png.
      approval: write
      output: ./schemas/time-roc-result.json
      timeoutSeconds: 600
---

# Clinical & Survival Biostatistics Workflow

Read [references/analysis-playbook.md](references/analysis-playbook.md) before choosing survival endpoints, covariate sets, or landmark horizons.

## Standard Workflow

```
clinical (+ optional omics matrix)
  ──▶ clin_inspect   : verify event encoding (1=event, 0=censored), time units, reverse-KM follow-up & EPV
  ──▶ clin_table1    : baseline cohort table (normality-routed continuous + Fisher/chi2 categorical + SMD)
  ──▶ clin_km        : Kaplan-Meier curves + Greenwood log-log CIs + exact log-rank + Number-at-Risk subplot
  ──▶ clin_cox       : univariable screen + multivariable Breslow Cox + Harrell's C-index + Schoenfeld PH test + forest plot
  ──▶ clin_time_roc  : out-of-fold L2-Cox risk score + IPCW time-dependent ROC at clinical landmarks (e.g. 1/3/5-yr)
```

## Hard Red Lines

1. **Verify event polarity first**: confirm `1` = event (death/progression) and `0` = right-censored via `clin_inspect`. Reversing `0/1` inverts every HR and KM curve.
2. **No minimum-p cutpoint fishing**: use `median`, `tertile`, `quartile`, or pre-defined clinical cutoffs in `clin_km`.
3. **EPV >= 10 for unpenalized multivariable Cox**: warn whenever `events_per_parameter < 10`; switch to penalized `clin_time_roc` or fewer covariates when events are scarce.
4. **Check proportional hazards (`ph_pvalue`)**: report Schoenfeld residual test p-values alongside multivariable Cox HRs.
5. **Out-of-fold evaluation for multi-feature risk scores**: never report apparent in-sample AUC on training samples as validation performance.
