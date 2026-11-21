Purpose: find the regions most specifically accessible in each cluster, and
optionally contrast two clusters directly. Principle: marker_regions tests every
cluster against the rest (cluster-vs-rest); diff_test compares two clusters
peak-by-peak with a logistic regression giving log fold change and significance
(the scATAC analogue of differential expression). Both run on the peak matrix
when peaks were called, otherwise on tiles. Interpretation: read the per-cluster
marker tables as regulatory hypotheses; a region accessible in one cluster but
not others suggests a cell-type-specific enhancer/promoter. Pitfalls: the
p-values are DESCRIPTIVE — clusters were defined from the same data (double
dipping), treat the ranking as evidence not significance; diff_test takes cell
INDICES not labels, and pairwise diff is O(clusters^2) so it is opt-in (run_diff);
on random/synthetic data no region passes the thresholds and that is an honest
result, not a bug. marker_regions is always cluster-vs-rest even when run_diff
is off.