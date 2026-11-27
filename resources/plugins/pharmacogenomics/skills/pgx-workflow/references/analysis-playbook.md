# Pharmacogenomics Analysis Playbook

## 1. 4-Parameter Log-Logistic (Hill) Dose-Response Model (`pgx_dose_response`)

Given positive drug concentrations $d > 0$ and normalized fractional cell viability $y \in [0, 1.5]$:

$$y(d) = E_{\infty} + \frac{E_0 - E_{\infty}}{1 + (d / \text{EC}_{50})^h}$$

- $E_0$: upper asymptote (vehicle-level viability, typically $\approx 1.0$).
- $E_{\infty}$: lower asymptote (residual viability at saturating concentration; maximum inhibition is $E_{\max} = 1 - E_{\infty}$).
- $\text{EC}_{50}$ (relative $\text{IC}_{50}$): inflection concentration halfway between $E_0$ and $E_{\infty}$.
- $h > 0$: Hill slope measuring cooperativity/steepness.
- **Absolute $\text{IC}_{50}$**: concentration where predicted viability crosses $y = 0.50$:
  $$\text{IC}_{50} = \text{EC}_{50} \left(\frac{E_0 - 0.50}{0.50 - E_{\infty}}\right)^{1/h}$$
  Defined only when $E_0 > 0.50 > E_{\infty}$. If $E_{\infty} \ge 0.50$ or $\text{IC}_{50} > \max(d)$, `beyond_tested_range` is set to `true` and `ic50_censored` is capped at `max_dose`.
- **Normalized AUC & AAC**: trapezoidal area under $y(\log_{10} d)$ divided by $(\log_{10} d_{\max} - \log_{10} d_{\min})$, clipped to $[0, 1]$, with $\text{AAC} = 1 - \text{AUC}$. Unlike $\text{IC}_{50}$, `AUC`/`AAC` remain well-defined even for resistant cell lines where 50% killing is never achieved.

## 2. Drug Combination Synergy Reference Models (`pgx_synergy`)

Convert fractional viability $y \in [0, 1]$ to fractional inhibition $E = 1 - y \in [0, 1]$ (reported on percent scale $\times 100$):

- **HSA (Highest Single Agent)**: $E_{\text{HSA}}(d_A, d_B) = \max(E_A(d_A), E_B(d_B))$.
- **Bliss Independence**: $E_{\text{Bliss}}(d_A, d_B) = E_A(d_A) + E_B(d_B) - E_A(d_A) E_B(d_B)$.
- **Loewe Additivity**: expected inhibition $E_{\text{Loewe}}$ solving the isobologram equation $\frac{d_A}{D_A(E)} + \frac{d_B}{D_B(E)} = 1$ using the fitted monotherapy Hill curves.
- **ZIP (Zero Interaction Potency)**: fits conditional Hill curves of Drug A in the presence of each fixed $d_B$ (and Drug B at each fixed $d_A$) and compares smoothed combined potency against $E_{\text{Bliss}}^{\text{fitted}}$.
- Synergy score $\delta = 100 \times (E_{\text{observed}} - E_{\text{expected}})$:
  - $\delta > +10$: **Synergistic**
  - $-10 \le \delta \le +10$: **Additive**
  - $\delta < -10$: **Antagonistic**

## 3. Biomarker Polarity Conventions (`pgx_biomarker`)

`pgx_biomarker` reports standardized `sensitivity_effect` where **positive values always indicate increased drug sensitivity** and **negative values indicate drug resistance**:
- When `higher_is_sensitive = false` (default for `IC50`, `LN_IC50`, `AUC`): `sensitivity_effect = -(mean_altered - mean_wt)` (or `-spearman_rho`).
- When `higher_is_sensitive = true` (for `AAC`, `% inhibition`): `sensitivity_effect = +(mean_altered - mean_wt)` (or `+spearman_rho`).

## 4. High-Throughput Screen Plate QC (`pgx_hts_qc`)

- **Zhang's $Z'$-factor**: $Z' = 1 - \frac{3(\sigma_+ + \sigma_-)}{|\mu_+ - \mu_-|}$. Plates with $Z' \ge 0.5$ are **excellent**, $0 \le Z' < 0.5$ are **marginal**, and $Z' < 0$ are **failed** and should be excluded or repeated before secondary dose-response follow-up.
- **Tukey Median Polish `B-score`**: fits additive row ($\alpha_r$) and column ($\beta_c$) medians on sample wells (masking positive/negative controls so extreme control values do not distort row/column estimates), then normalizes residuals by $1.4826 \cdot \text{MAD}$. Enable `--apply_bscore` whenever `row_anova_p < 0.01`, `col_anova_p < 0.01`, or `edge_vs_interior_p < 0.01`.

## 5. Cross-Cohort Drug Sensitivity Imputation (`pgx_predict_cohort`)

- Harmonizes gene-wise mean and variance between preclinical cell lines and patient tumors (applying empirical-Bayes shrinkage when clinical $n < 5$), then trains 5-fold cross-validated `RidgeCV` models per drug.
- Always check `drug_model_cv_metrics.tsv`: drugs with `cv_pearson_r < min_cv_r` carry `low_predictability_warning = true` and should not be prioritized as clinical stratification biomarkers without external validation.
