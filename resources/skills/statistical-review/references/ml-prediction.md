# Machine Learning & Clinical Prediction Audit Branch

Load this reference ONLY when auditing machine learning, deep learning, biomarker signatures, or clinical risk prediction models.

## 1. Domain Redlines Matrix

| Flaw | Why It Breaks Inference | Required Fix |
| --- | --- | --- |
| **Preprocessing / Feature Leakage** | Scaling, imputation, PCA, SMOTE, or feature selection on the full dataset before train/test split inflates `AUROC` | Fit all preprocessing, feature selection, and resampling strictly **inside each training fold** (`Pipeline` / nested CV) |
| **Patient / Donor Split Leakage** | Splitting multiple slides, biopsies, or cells from the same patient across train and test sets | Enforce **GroupKFold / patient-level split** (`donorBySample`) |
| **Discrimination Without Calibration** | High `AUROC` can still produce severely miscalibrated risk probabilities | Report **calibration curve / slope / intercept / Brier score** alongside `AUROC` |
| **Imbalanced-Class Metric Misuse** | Reporting only `Accuracy` or `ROC-AUC` on 2% positive prevalence hides poor precision | Report **PR-AUC**, sensitivity/specificity at pre-specified threshold, and **Decision Curve Analysis (DCA)** |
| **No External Validation** | Internal CV alone cannot measure domain shift across centers/scanners/batches | Validate on an **independent external cohort** (or temporal/center split) per **TRIPOD+AI** |

## 2. Worked Audit Micro-Example

- *Flawed Text*: *"We selected the top 15 genes by differential expression across all 120 patients, balanced classes with SMOTE, and achieved 5-fold CV `AUROC = 0.94`."*
- *Audit Findings*:
  1. **Critical — Dual data leakage**: Both feature selection and SMOTE occurred *outside* CV folds, leaking test-fold labels into training. Move DE feature selection and resampling strictly inside each training fold (**nested CV**).
  2. **Major — Missing calibration & external validation**: Add Brier score, calibration plot, `95% CI` via bootstrap, and an independent external cohort (or clearly label as internal discovery only).
