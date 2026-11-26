---
name: prediction-causal-workflow
description: >-
  Clinical prediction validation, diagnostic discrimination, causal inference,
  nonlinear dose-response, and evidence synthesis workflow in Python: DeLong
  ROC & paired AUC comparison + PR-AUC + Youden/fixed-cutoff 2x2 metrics;
  TRIPOD calibration (Brier, IPA, CITL, slope, HL) & Vickers Decision Curve
  Analysis (DCA); logistic propensity score positivity audit, 1:k caliper PSM,
  stabilized IPTW with ESS, Love-plot SMD audit & ATE/ATT estimation; Harrell
  Restricted Cubic Splines (RCS) with overall/non-linearity Wald tests; and
  Fixed/DL/PM/HKSJ meta-analysis with Egger's test, subgroups, and LOO. Use for
  诊断ROC, DeLong检验, 校准曲线, 临床决策曲线, DCA, 倾向性评分, PSM, IPTW,
  限制性立方样条, RCS, Meta分析, 森林图, or 漏斗图.
phi:
  environment: plugin:clinical
  attachTo:
    - Clinical
  scripts:
    - name: diagnostic_roc
      description: >-
        Evaluate 1 or more continuous diagnostic biomarkers / risk scores against
        a binary clinical outcome: empirical ROC AUC with exact Sun & Xu / DeLong
        (1988) 95% CIs, paired two-sided DeLong Z-tests with BH FDR comparing each
        predictor to the reference predictor, PR-AUC with bootstrap 95% CI, and
        2x2 sensitivity/specificity/PPV/NPV/LR+/LR- with Wilson 95% CIs.
      run: [python, ./scripts/clin_diagnostic_roc.py]
      args:
        type: object
        additionalProperties: false
        required: [clinical, outcome_col, predictors, output_dir]
        properties:
          clinical:
            type: string
            format: input-path
            description: Clinical dataset TSV/CSV.
          outcome_col:
            type: string
            description: Binary outcome column (e.g. 0/1 or Control/Case).
          predictors:
            type: array
            items:
              type: string
            description: Continuous biomarker or predicted risk score columns (first is reference for paired DeLong).
          positive_label:
            type: string
            description: Optional explicit value in outcome_col representing the positive class.
          threshold:
            type: number
            description: Optional pre-specified clinical cutoff (if omitted, uses Youden's J with resubstitution caveat).
          n_bootstrap:
            type: integer
            minimum: 50
            maximum: 5000
            description: Stratified bootstrap replicates for PR-AUC 95% CI (default 500).
          seed:
            type: integer
            description: Random seed for bootstrap (default 42).
          output_dir:
            type: string
            format: project-path
            description: Output directory for diagnostic_metrics.tsv, delong_comparisons.tsv, and diagnostic_roc_pr.png.
      approval: write
      output: ./schemas/diagnostic-roc-result.json
      timeoutSeconds: 300
    - name: calibration_dca
      description: >-
        TRIPOD calibration diagnostics and Vickers Decision Curve Analysis (DCA)
        for predicted probabilities (prob_col) or leak-free stratified K-fold
        out-of-fold L2 logistic regression on features: computes Brier score,
        scaled Brier (IPA), E/O ratio, calibration-in-the-large, calibration
        slope, Hosmer-Lemeshow test, and DCA net benefit / interventions avoided.
      run: [python, ./scripts/clin_calibration_dca.py]
      args:
        type: object
        additionalProperties: false
        required: [clinical, outcome_col, output_dir]
        properties:
          clinical:
            type: string
            format: input-path
            description: Clinical dataset TSV/CSV.
          outcome_col:
            type: string
            description: Binary outcome column.
          prob_col:
            type: string
            description: Pre-computed predicted probability column in [0, 1] (omit when passing features).
          features:
            type: array
            items:
              type: string
            description: Predictor columns for leak-free stratified K-fold out-of-fold logistic regression.
          positive_label:
            type: string
            description: Optional explicit positive outcome label.
          n_bins:
            type: integer
            minimum: 4
            maximum: 20
            description: Number of calibration quantile bins (default 10).
          n_folds:
            type: integer
            minimum: 3
            maximum: 10
            description: Cross-validation folds when fitting features (default 5).
          seed:
            type: integer
            description: Random seed for stratified CV (default 42).
          output_dir:
            type: string
            format: project-path
            description: Output directory for calibration_bins.tsv, dca_curve.tsv, and calibration_dca.png.
      approval: write
      output: ./schemas/calibration-dca-result.json
      timeoutSeconds: 300
    - name: psm_iptw
      description: >-
        Propensity score causal inference via greedy 1:k caliper matching (PSM)
        and stabilized truncated Inverse Probability of Treatment Weighting
        (IPTW): checks positivity/common support, computes unweighted/PSM/IPTW
        covariate SMDs (flagging |SMD| >= 0.10), exports matched cohort, and
        estimates causal treatment effects for continuous, binary, or survival outcomes.
      run: [python, ./scripts/clin_psm_iptw.py]
      args:
        type: object
        additionalProperties: false
        required: [clinical, treatment_col, covariates, output_dir]
        properties:
          clinical:
            type: string
            format: input-path
            description: Clinical cohort table TSV/CSV.
          treatment_col:
            type: string
            description: Binary treatment/exposure column (0/1 or Control/Treated).
          covariates:
            type: array
            items:
              type: string
            description: Baseline confounder columns for propensity score modeling.
          treated_label:
            type: string
            description: Optional explicit value in treatment_col denoting the treated group.
          outcome_col:
            type: string
            description: Optional outcome column (or survival event column) for causal effect estimation.
          outcome_type:
            type: string
            enum: [auto, binary, continuous, survival]
            description: Outcome type (auto-inferred by default).
          time_col:
            type: string
            description: Follow-up time column (required when outcome_type is survival).
          method:
            type: string
            enum: [both, psm, iptw]
            description: Propensity score adjustment method (default both).
          caliper_sd:
            type: number
            description: Matching caliper in units of SD(logit PS) (default 0.2).
          ratio:
            type: integer
            minimum: 1
            maximum: 5
            description: Maximum control matches per treated subject in 1:k PSM (default 1).
          truncate_pct:
            type: number
            description: Lower/upper quantile truncation for stabilized IPTW weights (default 0.01).
          sample_col:
            type: string
            description: Optional sample/patient ID column.
          output_dir:
            type: string
            format: project-path
            description: Output directory for covariate_balance_smd.tsv, psm_matched_cohort.tsv, and psm_iptw_diagnostics.png.
      approval: write
      output: ./schemas/psm-iptw-result.json
      timeoutSeconds: 600
    - name: rcs
      description: >-
        Fit Harrell's Restricted Cubic Spline (RCS) regression (3, 4, or 5 knots)
        for linear, logistic, or Cox PH models with optional covariate adjustment:
        performs Overall Association and Non-Linearity Wald chi-square tests and
        computes delta-method 95% CI dose-response curves relative to ref_value.
      run: [python, ./scripts/clin_rcs.py]
      args:
        type: object
        additionalProperties: false
        required: [clinical, exposure_col, outcome_col, output_dir]
        properties:
          clinical:
            type: string
            format: input-path
            description: Clinical dataset TSV/CSV.
          exposure_col:
            type: string
            description: Continuous exposure column to model via restricted cubic splines.
          outcome_col:
            type: string
            description: Outcome column (or survival event indicator for Cox).
          model_type:
            type: string
            enum: [auto, logistic, cox, linear]
            description: Regression family (auto-inferred by default; cox requires time_col).
          time_col:
            type: string
            description: Follow-up time column (required for Cox RCS).
          covariates:
            type: array
            items:
              type: string
            description: Optional confounding covariates to adjust for.
          n_knots:
            type: integer
            enum: [3, 4, 5]
            description: Number of Harrell canonical spline knots (default 4).
          ref_value:
            type: number
            description: Reference exposure value where OR/HR=1 or delta=0 (default cohort median).
          output_dir:
            type: string
            format: project-path
            description: Output directory for rcs_curve.tsv and rcs_dose_response.png.
      approval: write
      output: ./schemas/rcs-result.json
      timeoutSeconds: 300
    - name: meta
      description: >-
        Study-level Fixed-Effect and Random-Effects (DerSimonian-Laird,
        Paule-Mandel, and Hartung-Knapp-Sidik-Jonkman) meta-analysis from
        pre-computed effect sizes/SEs/CIs or 2x2 contingency counts (with
        Haldane-Anscombe 0.5 continuity correction): reports Q, I^2, tau^2, 95%
        prediction interval, Egger's test, subgroup Q_between, Leave-One-Out
        sensitivity, and forest + funnel plots.
      run: [python, ./scripts/clin_meta.py]
      args:
        type: object
        additionalProperties: false
        required: [studies, study_col, output_dir]
        properties:
          studies:
            type: string
            format: input-path
            description: Study-level summary table TSV/CSV.
          study_col:
            type: string
            description: Study identifier column.
          effect_col:
            type: string
            description: Study effect estimate column (log-HR/log-OR/SMD/MD, or ratio-scale HR/OR/RR with is_ratio_scale).
          se_col:
            type: string
            description: Standard error column (on log scale when is_ratio_scale is used).
          ci_low_col:
            type: string
            description: Lower 95% CI column (alternative to se_col).
          ci_high_col:
            type: string
            description: Upper 95% CI column (alternative to se_col).
          is_ratio_scale:
            type: boolean
            description: Set true when effect_col and ci_low/ci_high are exponentiated HR/OR/RR > 0 (default false).
          ai_col:
            type: string
            description: 2x2 treated events count column.
          bi_col:
            type: string
            description: 2x2 treated non-events count column.
          ci_col:
            type: string
            description: 2x2 control events count column.
          di_col:
            type: string
            description: 2x2 control non-events count column.
          measure:
            type: string
            enum: [OR, RR, RD]
            description: Effect measure when pooling 2x2 counts (default OR).
          subgroup_col:
            type: string
            description: Optional study-level subgroup column for stratified pooling and Q_between test.
          output_dir:
            type: string
            format: project-path
            description: Output directory for meta_summary.tsv, meta_leave_one_out.tsv, and meta_forest_funnel.png.
      approval: write
      output: ./schemas/meta-result.json
      timeoutSeconds: 300
---

# Prediction Validation, Causal Inference & Meta-Analysis Workflow

Read [references/analysis-playbook.md](references/analysis-playbook.md) before interpreting DeLong comparisons, TRIPOD calibration/DCA, PSM/IPTW balance, RCS non-linearity, or meta-analysis heterogeneity.

## Standard Workflows

```
Diagnostic / Prediction Validation:
  ──▶ clin_diagnostic_roc   : DeLong AUC + paired DeLong Z-test + PR-AUC + Youden/fixed 2x2 Wilson CIs
  ──▶ clin_calibration_dca  : TRIPOD calibration (Brier, IPA, CITL, slope, HL) + Vickers DCA net benefit

Observational Causal Inference & Dose-Response:
  ──▶ clin_psm_iptw         : positivity audit + 1:k caliper PSM + stabilized IPTW + |SMD| < 0.1 + ATE/ATT
  ──▶ clin_rcs              : Harrell 3/4/5-knot Restricted Cubic Splines + overall & non-linearity Wald tests

Evidence Synthesis:
  ──▶ clin_meta             : Fixed + DL/PM/HKSJ random-effects pooling + I^2/tau^2 + 95% PI + Egger + LOO
```

## Hard Red Lines

1. **Report calibration and clinical utility alongside AUC**: high ROC AUC alone does not prove clinical utility; always run `clin_calibration_dca` for TRIPOD calibration slope, calibration-in-the-large, and DCA net benefit.
2. **Never report apparent Youden's $J$ metrics as validated cutoffs**: unless `threshold` was externally pre-specified, label `clin_diagnostic_roc` sensitivity/specificity as resubstitution estimates requiring external validation.
3. **Check positivity and `|SMD| < 0.10` before causal claims**: in `clin_psm_iptw`, flag extreme propensity scores (`<0.05` or `>0.95`), low ESS, or residual `|SMD| >= 0.10`.
4. **Do not assume log-linearity without testing**: run `clin_rcs` and report `p_nonlinear` before modeling continuous exposures linearly; use 4 knots unless sample size is small (`n < 100` -> 3 knots).
5. **Respect heterogeneity & publication-bias power in meta-analysis**: report the 95% prediction interval alongside the pooled RE CI when $I^2$ is substantial, and warn that Egger's test is underpowered when $k < 10$.
