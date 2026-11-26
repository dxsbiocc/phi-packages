# Clinical & Survival Biostatistics Playbook

Decision manual for the `Clinical` agent (`clin_*` survival tools).

## Stage 1 — Cohort & Follow-up Audit (`clin_inspect`)
- **Call**: `clin_inspect(clinical=..., time_col=..., event_col=..., matrix=...)`
- **Read**: `n_valid_survival`, `n_events`, `censoring_fraction`, `inferred_time_unit`, `median_followup_reverse_km`, `epv_max_parameters_at_10`.
- **Decide**:
  - Choose `landmarks` strictly below `max_event_time` where sufficient uncensored patients remain at risk (e.g., `12, 36, 60` when `inferred_time_unit == "months"`, or `1, 3, 5` when `"years"`).
  - Cap unpenalized multivariable Cox degrees of freedom at `epv_max_parameters_at_10` (`floor(n_events / 10)`).
- **Red line**: Drop or investigate records with `time <= 0` (immortal-time / clerical error) rather than silently recoding them to positive numbers.

## Stage 2 — Baseline Characteristics (`clin_table1`)
- **Call**: `clin_table1(clinical=..., variables=[...], group_col=..., output_dir=...)`
- **Read**: `variables` (`test`, `pvalue`, `smd`, `summary_by_group`).
- **Decide**: Variables with baseline imbalance (`pvalue < 0.10` or `|smd| > 0.20` between exposure groups) are primary confounding candidates to adjust for in `clin_cox` or `clin_psm_iptw`.

## Stage 3 — Kaplan-Meier & Log-rank (`clin_km`)
- **Call**: `clin_km(clinical=..., time_col=..., event_col=..., group_col=... | biomarker_col=..., cutoff="median", output_dir=...)`
- **Read**: `logrank_pvalue`, `groups` (`n`, `events`, `median_survival`, `median_ci_low`, `median_ci_high`), `landmark_survival`.
- **Red line**: Crossing KM curves warn of non-proportional hazards — verify with Schoenfeld residuals in `clin_cox`. Never sweep arbitrary percentiles to minimize `logrank_pvalue`.

## Stage 4 — Cox Proportional Hazards (`clin_cox`)
- **Call**: `clin_cox(clinical=..., time_col=..., event_col=..., covariates=[...], mode="both", output_dir=...)`
- **Read**: `multivariable.c_index`, `multivariable.events_per_parameter`, `multivariable.terms` (`hr`, `ci_low`, `ci_high`, `pvalue`, `ph_pvalue`), `global_ph_violation`.
- **Red line**: If `events_per_parameter < 10`, flag overfitting risk explicitly. If any term has `ph_pvalue < 0.05`, state that the PH assumption is violated for that variable.

## Stage 5 — Time-Dependent ROC & Risk Signatures (`clin_time_roc`)
- **Call**: `clin_time_roc(clinical=..., time_col=..., event_col=..., score_col=... | feature_cols=[...], landmarks=[...], output_dir=...)`
- **Read**: `evaluation_mode` (`precomputed_score` vs `oof_cv_l2_cox`), `c_index`, `landmarks` (`horizon`, `n_cases`, `n_controls`, `auc`), `oof_km_logrank_pvalue`.
- **Red line**: Only landmarks with both `n_cases >= 3` (`T_i <= t` and `delta_i == 1`) and `n_controls >= 3` (`T_i > t`) yield reliable IPCW AUC estimates; skipped horizons are reported transparently.

## Stage 6 — Competing Risks & Fine-Gray (`clin_competing_risk`)
- **Call**: `clin_competing_risk(clinical=..., time_col=..., status_col=..., primary_code="1", group_col=..., covariates=[...], output_dir=...)`
- **Read**: `n_primary_events`, `n_competing_events`, `naive_1_minus_km_overestimation`, `grays_test` (`chisq`, `pvalue`), `fine_gray_model.terms` (`shr`, `ci_low`, `ci_high`, `pvalue`, `padj`).
- **Decide**:
  - Use **Aalen-Johansen CIF** and **Fine-Gray `sHR`** when modeling absolute cumulative risk / prognosis in the presence of competing events (e.g., disease-specific death vs non-disease death, relapse vs transplant).
  - Note the distinction: cause-specific Cox HR reflects instantaneous rate among event-free survivors (etiologic), whereas Fine-Gray `sHR` reflects cumulative incidence subdistribution effects (prognostic).
- **Red line**: Never interpret `1 - KM` as absolute event probability when `n_competing_events > 0`; report the exact `naive_1_minus_km_overestimation` bias.
