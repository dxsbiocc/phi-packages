# `mg_multivariable` — CLR multivariable differential abundance

## Purpose

Performs compositionally-aware, covariate-adjusted differential abundance
testing on any `features × samples` matrix (taxa or functional groups) by
applying a per-sample Centered Log-Ratio (CLR) transformation followed by OLS
linear regression (`CLR(feature) ~ exposure + covariates`) with exact Student-t
p-values and Benjamini-Hochberg FDR control.

## When to use

- When clinical or environmental confounders (e.g., `sex`, `age`, `bmi`,
  `batch`) must be adjusted alongside the primary `exposure` variable.
- When relative abundances suffer from compositional closure bias that can
  distort unadjusted Mann-Whitney U fold changes.

## Caveats

- Requires enough samples for residual degrees of freedom (`n > p + 1`).
- Fits a single Gaussian OLS model on CLR-transformed values; if extreme zero
  inflation requires a two-part hurdle model, use R-based MaAsLin 3 or ANCOM-BC.
