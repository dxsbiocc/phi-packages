# Prediction Validation, Causal Inference & Meta-Analysis Playbook

Decision manual for the `Clinical` agent (`prediction-causal-workflow` tools).

## 1. Diagnostic ROC & Paired DeLong Comparison (`clin_diagnostic_roc`)
- **TRIPOD / STARD alignment**:
  - Report empirical ROC AUC with exact **Sun & Xu / DeLong (1988)** structural-components 95% CI, plus **Precision-Recall AUC (AP)** when outcome prevalence is imbalanced (`< 20%`).
  - Use **paired two-sided DeLong $Z$-tests** (`delong_comparisons`) when comparing candidate biomarkers or models on the exact same patients rather than comparing overlapping confidence intervals.
  - Inspect `direction` (`positive` vs `negative`): if a biomarker is protective (`raw_auc < 0.5`), the tool aligns polarity for discrimination comparison and flags `direction_inverted: true`.
- **Red line**: Selecting an optimal cutoff via Youden's $J$ on the evaluation set introduces optimism (resubstitution bias). Report it as exploratory unless `threshold` was pre-specified externally.

## 2. TRIPOD Calibration & Decision Curve Analysis (`clin_calibration_dca`)
- **Leak-free probability evaluation**:
  - Pass `prob_col` when validating an external/locked model, or pass `features` to generate leak-free **stratified K-fold out-of-fold** predicted probabilities.
- **Calibration hierarchy (Van Calster et al. / TRIPOD)**:
  - **Calibration-in-the-large (CITL)**: ideal `0.0` (`< 0` = over-prediction, `> 0` = under-prediction).
  - **Calibration slope**: ideal `1.0` (`< 1.0` = extreme/overfit predictions; `> 1.0` = under-extreme).
  - **Brier score & IPA (Scaled Brier)**: `IPA = 1 - Brier / Brier_null` (`> 0` beats prevalence-only model).
  - **Hosmer-Lemeshow test**: sensitive to large $N$ — interpret alongside the visual calibration plot and $E/O$ ratio.
- **Vickers & Elkin (2006) Decision Curve Analysis (DCA)**:
  - Net Benefit $\text{NB}(p_t) = \frac{\text{TP}}{N} - \frac{\text{FP}}{N}\frac{p_t}{1 - p_t}$ measures true positives gained net of false positives weighted by threshold odds.
  - A model is clinically useful only across `useful_threshold_range` where $\text{NB}_{\text{model}}(p_t) > \max(\text{NB}_{\text{Treat-All}}(p_t), 0)$.

## 3. Propensity Score Matching & IPTW (`clin_psm_iptw`)
- **Rosenbaum-Rubin assumptions**:
  1. **Positivity / common support**: inspect `positivity_audit` (`n_extreme_ps_lt05`, `n_extreme_ps_gt95`, `overlap_range`). Severe non-overlap breaks causal identifiability.
  2. **Covariate balance (`|SMD| < 0.10`)**: inspect `covariate_balance_smd.tsv` and the Love plot. Any variable with post-adjustment `|SMD| >= 0.10` represents residual confounding.
  3. **Effective Sample Size (ESS)**: large IPTW weights shrink ESS (`ess_treated`, `ess_control`); stabilized + truncated weights (`truncate_pct=0.01`) control variance inflation.
- **Estimand distinction**:
  - **PSM** estimates the **ATT** (Average Treatment Effect on the Treated) in matched pairs within `caliper_sd * SD(logit PS)`.
  - **Stabilized IPTW** estimates the **ATE** (Average Treatment Effect across the population) with Huber-White sandwich 95% CIs.

## 4. Harrell's Restricted Cubic Splines (`clin_rcs`)
- **Knot placement (Harrell *Regression Modeling Strategies*)**:
  - `3` knots: `10, 50, 90` percentiles (use when `n < 100` or `events < 40`).
  - `4` knots (default): `5, 35, 65, 95` percentiles.
  - `5` knots: `5, 27.5, 50, 72.5, 95` percentiles (use for large cohorts `n >= 500`).
- **Interpretation**:
  - Report both **`p_overall`** ($H_0: \text{all spline terms} = 0$, $df = K - 1$) and **`p_nonlinear`** ($H_0: \text{nonlinear terms} = 0$, $df = K - 2$).
  - Only claim a U-shaped, J-shaped, or threshold dose-response when `p_nonlinear < 0.05`.

## 5. Meta-Analysis (`clin_meta`)
- **Heterogeneity & Pooling**:
  - Report **Fixed-Effect** alongside **Random-Effects (DerSimonian-Laird / Paule-Mandel)** and **Hartung-Knapp-Sidik-Jonkman (HKSJ)** 95% CIs (more robust when $k$ is small).
  - Always inspect **$I^2$**, **$\tau^2$**, and the **95% Prediction Interval** ($k \ge 3$): unlike the pooled CI (uncertainty of the mean effect), the prediction interval shows the range of true effects expected in a future study.
- **Bias & Sensitivity**:
  - Check **Leave-One-Out** (`meta_leave_one_out.tsv`) to verify no single study drives significance.
  - **Egger's regression asymmetry test** requires $k \ge 3$ and has low statistical power when $k < 10$.
