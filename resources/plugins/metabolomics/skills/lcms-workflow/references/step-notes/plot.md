# mb_plot — clustered feature heatmap

## Principle

A heatmap of the top-variance features shows sample structure at a glance:
groups, outliers, batch stripes. Deterministic rendering; the agent chooses
what to display.

## Choices

- `top` (default 50): features ranked by variance across samples. Variance
  on raw intensities mostly ranks abundance — log first when that is not
  the intent.
- `trans log10` + `pseudocount` (default 1e-6): compresses dynamic range;
  missing values are filled with 0 before transform (pseudocount keeps it
  finite). Negative values are refused.
- `cluster`: average-linkage hierarchical clustering (Euclidean on the
  displayed values) of rows, columns, or both. Clustering raw intensities
  groups by abundance level; pair with log10 for pattern-based grouping.

## Interpretation

- Column clusters that disagree with the metadata groups are a finding:
  batch effects, mislabeled samples, or real biology — investigate with
  PCA (mb_multivariate) before claiming either.
- Row blocks are co-varying feature sets; they are not pathways. Annotate
  (mb_annotate) and enrich (mb_pathway) before functional claims.
- The color scale is per-figure; do not compare colors across figures.

## Pitfalls

- Fewer than 2 sample columns → the tool fails; a heatmap of one sample is
  meaningless.
- All-missing features are dropped; the count shown is `n_features_shown`.
- Score plots and volcanoes come from mb_multivariate / mb_univariate, not
  here. Every figure carries a `.phi-artifact.json` descriptor.