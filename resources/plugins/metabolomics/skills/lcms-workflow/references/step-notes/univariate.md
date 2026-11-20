# mb_univariate — per-feature tests, FDR, fold change

## Principle

One test per feature, then multiplicity control across features. The test
is a choice: `wilcoxon` (Mann-Whitney U, two-sided, tie-corrected) is the
robust default for skewed intensities; `welch` (unequal-variance t) has
more power on roughly normal data (e.g. after log); `kruskal` handles >2
groups without a case/control contrast.

## Fold change

- `logged: true`: log2FC = mean(case) − mean(control), FC = 2^log2FC.
- Raw intensities: FC = mean(case)/mean(control); zero control mean with
  positive case mean → inf (reported, not hidden); both zero → FC 1.
- Never quote a raw ratio as if it were a log difference.

## Multiplicity and thresholds

- BH FDR across all tested features (monotone step-up on sorted p).
- Significance = FDR < threshold AND |log2FC| ≥ fc_threshold (two-group);
  FDR only for kruskal. Defaults 0.05 / 1 are conventions — pass explicitly.
- Features with <2 observed values in a group are skipped and counted in
  `n_tested` vs `n_features`.

## Interpretation

- With n=3/group the smallest possible two-sided Wilcoxon p is 0.1: no
  feature can pass FDR 0.05. A null result is then about power, not biology.
- Volcano: x = log2FC, y = −log10 p; red = passes both thresholds. The
  dashed lines are the thresholds you passed, not universal cutoffs.
- Intersect with VIP > 1 (mb_multivariate) for the conventional
  "multivariate + univariate" candidate list.

## Pitfalls

- Run on biological samples only; QC columns in the matrix will fail
  metadata coverage or pollute groups.
- Constant features (all equal in both groups) get p = 1 by definition,
  not a scipy NaN.