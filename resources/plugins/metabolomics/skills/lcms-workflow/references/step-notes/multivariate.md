# mb_multivariate — PCA, PLS-DA, VIP, permutation Q2

## Principle

PCA (SVD on centered/scaled data) is unsupervised: it shows the dominant
variance structure — batches, outliers, gradients — without labels.
PLS-DA (NIPALS PLS2 against one-hot group Y) is supervised: it finds
directions that maximize covariance with the labels. Supervised separation
is easy to fake; validation is the point (Westerhuis 2008).

## Scaling

- `uv`: center + divide by SD — every feature weighs equally; noise
  features get promoted.
- `pareto`: center + divide by sqrt(SD) — middle ground, common in
  metabolomics.
- `none`: centering only — abundant features dominate.
- Zero-variance features are centered to zero, never divided by zero.

## Model quality

- R2Y cum: fraction of label variance explained in-sample. High R2Y alone
  means nothing — PLS-DA overfits readily.
- Q2 cum: cross-validated predictive ability (LOO when n ≤ 50, seeded
  7-fold above). Scaling is recomputed inside each training fold.
- Permutation: labels shuffled `n_perm` times (seeded), full refit + CV
  each time; `q2_p` / `r2y_p` = (1 + #null ≥ observed) / (1 + n_perm).
  Q2 ≤ 0 or non-significant permutation p → overfit; report it.
- `n_components` is clamped to min(requested, n−2, p); the JSON says so.

## VIP

VIP_j weights each feature's contribution across components by the
explained Y sum of squares. Convention: VIP > 1 = important. Use it as a
ranking to intersect with univariate FDR/FC — never as sole evidence.

## Pitfalls

- Drop QC/blank columns before fitting; they are not biological groups.
- Small n makes Q2 noisy; a negative Q2 with n=6 is a power statement.
- OPLS-DA (ropls/mixOmics) is a gap. Do not relabel PLS-DA output as OPLS-DA.
- Same seed → identical permutations and folds; cite the seed in reports.