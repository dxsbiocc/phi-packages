# Survival Analysis & Clinical Epidemiology Audit Branch

Load this reference ONLY when auditing survival models, clinical cohorts, RCTs, or epidemiological studies.

## 1. Domain Redlines Matrix

| Flaw | Why It Breaks Inference | Required Fix |
| --- | --- | --- |
| **Unchecked Proportional Hazards (PH)** | Time-varying hazard makes a single `HR` misleading | Test **Schoenfeld residuals** (`cox.zph`); if violated, use stratified Cox, time-interaction, or **RMST** |
| **Ignoring Competing Risks** | Standard KM/Cox censoring death from other causes overestimates disease-specific incidence | **Aalen-Johansen** CIF + **Fine-Gray** subdistribution hazard (or cause-specific Cox with explicit framing) |
| **Immortal Time Bias** | Defining responder/treatment group *after* baseline `t = 0` forces responders to survive until classification | **Landmark analysis** or **time-dependent Cox** covariate |
| **Naive Follow-Up Time** | Reporting median survival time or raw mean follow-up underestimates true follow-up under censoring | Report **reverse Kaplan-Meier** median follow-up (`IQR`) |
| **Overfitting (`EPV < 10`)** | Fitting too many covariates per event destabilizes `HR`/`OR` | Pre-specify fewer confounders, use **Firth / ridge / LASSO penalization**, and run sensitivity checks |
| **Subgroup Fallacy** | Claiming subgroup A benefits because `P = 0.02` in A vs `P = 0.18` in B | Run a formal **interaction test** (`treatment * subgroup`) |
| **Ratio & Pooling Paradoxes** | Correlating `X/Z` vs `Y/Z` (spurious ratio correlation); pooling unadjusted cohorts (**Simpson's paradox**) | Regress `Y ~ X + Z` (ANCOVA/offset); stratify or adjust for center/cohort |

## 2. Worked Audit Micro-Example

- *Flawed Text*: *"Patients achieving complete response within 6 months had longer OS from diagnosis (`HR = 0.31`, `P = 0.001`, Cox model adjusted for 8 variables, `n = 48`, 22 deaths)."*
- *Audit Findings*:
  1. **Critical — Immortal time bias**: Response is defined post-baseline; switch to a **6-month landmark** or **time-dependent Cox** model.
  2. **Critical — Severe overfitting (`EPV = 22/8 = 2.75 < 10`)**: Limit to `<= 2` pre-specified confounders or use penalized Cox.
  3. **Major — Unreported PH & follow-up**: `[UNVERIFIABLE FROM TEXT — INSPECT CODE/RAW DATA]` — add Schoenfeld `P`-values and reverse-KM median follow-up.
