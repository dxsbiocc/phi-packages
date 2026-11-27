---
name: prediction-causal-workflow
description: >-
  Clinical prediction validation, diagnostic discrimination, causal inference,
  nonlinear dose-response, trial sample-size design, subgroup interaction,
  nomogram construction, and evidence synthesis workflow in Python: DeLong ROC &
  paired AUC comparison + PR-AUC + Youden/fixed-cutoff 2x2 metrics; TRIPOD
  calibration (Brier, IPA, CITL, slope, HL) & Vickers Decision Curve Analysis
  (DCA); logistic propensity score positivity audit, 1:k caliper PSM, stabilized
  IPTW with ESS, Love-plot SMD audit & ATE/ATT estimation; Harrell Restricted
  Cubic Splines (RCS); Fixed/DL/PM/HKSJ meta-analysis; Schoenfeld/Fleiss/Farrington-Manning/Hanley-McNeil
  sample size & power; formal multiplicative subgroup P-for-interaction forest
  plots; and Cox/Logistic 0-100 point clinical nomograms. Use for 诊断ROC,
  DeLong检验, 校准曲线, 临床决策曲线, DCA, 倾向性评分, PSM, IPTW, 限制性立方样条, RCS,
  Meta分析, 样本量计算, 检验效能, 亚组交互作用, 列线图, Nomogram, 森林图, or 漏斗图.
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
    - name: sample_size
      description: >-
        Prospective clinical study sample-size and statistical power calculator
        across survival_logrank (Schoenfeld 1983), two_proportions (Fleiss/Casagrande
        superiority and Farrington-Manning 1990 non-inferiority), two_means
        (exact non-central t-distribution), and diagnostic_auc (Hanley & McNeil
        1982), with dropout inflation and conservative/target/optimistic power curves.
      run: [python, ./scripts/clin_sample_size.py]
      args:
        type: object
        additionalProperties: false
        required: [design, output_dir]
        properties:
          design:
            type: string
            enum: [survival_logrank, two_proportions, two_means, diagnostic_auc]
            description: Clinical study endpoint / test family.
          output_dir:
            type: string
            format: project-path
            description: Output directory for power_curve_grid.tsv and sample_size_power_curve.png.
          solve_for:
            type: string
            enum: [sample_size, power]
            description: Solve for required sample size given target power, or achieved power given n_total (default sample_size).
          alpha:
            type: number
            description: Type I error rate alpha (default 0.05).
          power:
            type: number
            description: Target statistical power 1 - beta when solve_for is sample_size (default 0.80).
          n_total:
            type: integer
            minimum: 4
            description: Total enrolled sample size when solve_for is power.
          allocation_ratio:
            type: number
            description: Experimental-to-control arm allocation ratio k = N_exp / N_ctrl (default 1.0).
          one_sided:
            type: boolean
            description: Use one-sided significance testing instead of two-sided (default false).
          dropout_rate:
            type: number
            description: Anticipated attrition/dropout fraction in [0, 1) used to inflate enrolled N (default 0.10).
          hr:
            type: number
            description: Target hazard ratio (HR != 1) for survival_logrank.
          event_rate_ctrl:
            type: number
            description: Control-arm cumulative event probability in (0, 1] for survival_logrank.
          median_surv_ctrl:
            type: number
            description: Control-arm median survival time (alternative to event_rate_ctrl when paired with accrual_time/followup_time).
          accrual_time:
            type: number
            description: Uniform patient accrual duration in same time unit as median_surv_ctrl.
          followup_time:
            type: number
            description: Minimum post-accrual follow-up duration in same time unit as median_surv_ctrl.
          p_ctrl:
            type: number
            description: Control-arm event/response proportion in (0, 1) for two_proportions.
          p_exp:
            type: number
            description: Experimental-arm event/response proportion in (0, 1) for two_proportions.
          non_inferiority_margin:
            type: number
            description: Optional Farrington-Manning non-inferiority margin delta in (0, 1) for two_proportions.
          mean_diff:
            type: number
            description: Expected difference in means (mu_exp - mu_ctrl) for two_means.
          sd:
            type: number
            description: Common within-group standard deviation (> 0) for two_means.
          auc_target:
            type: number
            description: Expected diagnostic ROC AUC under H1 in (0, 1) for diagnostic_auc.
          auc_null:
            type: number
            description: Null hypothesis ROC AUC in (0, 1) for diagnostic_auc (default 0.50).
          prevalence:
            type: number
            description: Positive-case disease prevalence in (0, 1) for diagnostic_auc (default 0.30).
      approval: write
      output: ./schemas/sample-size-result.json
      timeoutSeconds: 300
    - name: subgroup_interaction
      description: >-
        Stratified subgroup effect estimation and formal multiplicative
        treatment-by-subgroup interaction testing (P for interaction via LRT /
        joint Wald chi-square) with BH-FDR multiplicity control and
        NEJM/Lancet-style subgroup forest plot across Cox PH, Logistic, and Linear models.
      run: [python, ./scripts/clin_subgroup_interaction.py]
      args:
        type: object
        additionalProperties: false
        required: [clinical, exposure_col, subgroups, outcome_col, output_dir]
        properties:
          clinical:
            type: string
            format: input-path
            description: Clinical cohort table TSV/CSV.
          exposure_col:
            type: string
            description: Primary treatment or exposure column (binary or continuous).
          subgroups:
            type: array
            items:
              type: string
            description: Candidate subgroup modifier columns (continuous columns are auto-dichotomized at the median).
          outcome_col:
            type: string
            description: Clinical outcome column (or survival event indicator when model_type is cox).
          model_type:
            type: string
            enum: [auto, cox, logistic, linear]
            description: Regression family (auto-inferred by default; cox requires time_col).
          time_col:
            type: string
            description: Follow-up time column (required for Cox PH subgroup interaction).
          adjust_covariates:
            type: array
            items:
              type: string
            description: Optional baseline covariates to adjust for across strata and interaction models.
          exposed_label:
            type: string
            description: Optional explicit value in exposure_col denoting the exposed/treated arm.
          output_dir:
            type: string
            format: project-path
            description: Output directory for subgroup_interaction.tsv and subgroup_forest_plot.png.
      approval: write
      output: ./schemas/subgroup-interaction-result.json
      timeoutSeconds: 300
    - name: nomogram
      description: >-
        Build a multivariable Cox PH (with Breslow baseline survival at landmark
        horizons) or Logistic clinical prediction alignment nomogram: exports a
        0-100 point scorecard per predictor, exact Total Points to predicted
        event probability / landmark survival lookup table, EPV check, and
        publication-ready nomogram figure.
      run: [python, ./scripts/clin_nomogram.py]
      args:
        type: object
        additionalProperties: false
        required: [clinical, predictors, outcome_col, output_dir]
        properties:
          clinical:
            type: string
            format: input-path
            description: Clinical dataset TSV/CSV.
          predictors:
            type: array
            items:
              type: string
            description: 2 to 12 clinical/molecular predictor columns included in the nomogram.
          outcome_col:
            type: string
            description: Binary outcome column (or survival event column when model_type is cox).
          model_type:
            type: string
            enum: [auto, cox, logistic]
            description: Nomogram regression family (auto selects cox when time_col is given, else logistic).
          time_col:
            type: string
            description: Follow-up time column (required for Cox survival nomogram).
          landmarks:
            type: array
            items:
              type: number
            description: Optional survival horizons (e.g. 12, 36, 60) mapped onto bottom survival scales in Cox mode.
          time_label:
            type: string
            description: Optional display unit label for landmark survival scales (e.g. months, years).
          output_dir:
            type: string
            format: project-path
            description: Output directory for nomogram_points_scorecard.tsv, nomogram_total_points_lookup.tsv, and clinical_nomogram.png.
      approval: write
      output: ./schemas/nomogram-result.json
      timeoutSeconds: 300
---

# Prediction Validation, Causal Inference, Trial Design & Meta-Analysis Workflow

Read [references/analysis-playbook.md](references/analysis-playbook.md) before interpreting DeLong comparisons, TRIPOD calibration/DCA, PSM/IPTW balance, RCS non-linearity, sample-size/power assumptions, subgroup interaction tests, nomograms, or meta-analysis heterogeneity.

## Standard Workflows

```
Study Design & Sample Size:
  ──▶ clin_sample_size          : Schoenfeld log-rank / Fleiss & Farrington-Manning proportions / exact nct / Hanley-McNeil AUC

Diagnostic / Prediction Modeling & Validation:
  ──▶ clin_nomogram             : Cox (landmark survival) or Logistic 0-100 point scorecard + alignment nomogram
  ──▶ clin_diagnostic_roc       : DeLong AUC + paired DeLong Z-test + PR-AUC + Youden/fixed 2x2 Wilson CIs
  ──▶ clin_calibration_dca      : TRIPOD calibration (Brier, IPA, CITL, slope, HL) + Vickers DCA net benefit

Observational Causal Inference, Subgroup Heterogeneity & Dose-Response:
  ──▶ clin_psm_iptw             : positivity audit + 1:k caliper PSM + stabilized IPTW + |SMD| < 0.1 + ATE/ATT
  ──▶ clin_subgroup_interaction : stratified forest plot + formal multiplicative P_interaction + BH-FDR
  ──▶ clin_rcs                  : Harrell 3/4/5-knot Restricted Cubic Splines + overall & non-linearity Wald tests

Evidence Synthesis:
  ──▶ clin_meta                 : Fixed + DL/PM/HKSJ random-effects pooling + I^2/tau^2 + 95% PI + Egger + LOO
```

## Hard Red Lines

1. **Report calibration and clinical utility alongside AUC / Nomograms**: high ROC AUC or C-index alone does not prove clinical utility; always pair `clin_nomogram` and `clin_diagnostic_roc` with `clin_calibration_dca` (calibration slope, CITL, and DCA net benefit) and verify `events_per_variable >= 10`.
2. **Never report apparent Youden's $J$ metrics as validated cutoffs**: unless `threshold` was externally pre-specified, label `clin_diagnostic_roc` sensitivity/specificity as resubstitution estimates requiring external validation.
3. **Check positivity and `|SMD| < 0.10` before causal claims**: in `clin_psm_iptw`, flag extreme propensity scores (`<0.05` or `>0.95`), low ESS, or residual `|SMD| >= 0.10`.
4. **Subgroup claims require formal $P_{\text{interaction}}$**: never claim treatment heterogeneity because one subgroup has $p < 0.05$ and another has $p > 0.05$; require `p_interaction < 0.05` (and `padj_interaction < 0.05` across multiple exploratory subgroups).
5. **Do not assume log-linearity without testing**: run `clin_rcs` and report `p_nonlinear` before modeling continuous exposures linearly; use 4 knots unless sample size is small (`n < 100` -> 3 knots).
6. **Respect heterogeneity & publication-bias power in meta-analysis**: report the 95% prediction interval alongside the pooled RE CI when $I^2$ is substantial, and warn that Egger's test is underpowered when $k < 10$.
