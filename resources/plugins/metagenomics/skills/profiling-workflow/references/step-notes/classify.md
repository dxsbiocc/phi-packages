# `mg_classify` — Leak-free cross-validated biomarker classification

## Purpose

Evaluates whether taxonomic (species/SGB) or functional (UniRef90 gene-family
or KO) profiles can discriminate a binary `case` vs `control` outcome using
stratified k-fold cross-validation with strictly in-fold log1p + z-score
scaling, in-fold feature selection (`mrmr` or `mannwhitney`), and L2-regularized
logistic regression.

## Outputs

- `oof-predictions.tsv`: out-of-fold predicted probabilities and Youden-optimal
  binary calls per sample.
- `feature-importance.tsv`: consensus features ranked by fold selection
  frequency and mean standardized logistic regression coefficient.
- `roc-curve.png`: pooled out-of-fold ROC curve with AUROC and Youden operating
  point (`sensitivity`, `specificity`).

## Interpretation and pitfalls

- **Zero leakage by design**: scaling and feature selection are fit solely on
  each training fold; never pre-filter features using full-cohort p-values
  before calling `mg_classify`, or AUROC will be optimistically biased.
- **Cohort size**: small cohorts (<30 samples per class) have wide fold-to-fold
  AUROC variance; always report `fold_auroc_sd` alongside pooled `oof_auroc`.
