# Composition figures

**Purpose.** Publication-starting views of the matrix: stacked barplots and
taxa heatmaps.

**Principle.** Top-N by mean abundance plus an "Other" bucket keeps the
legend readable; group-mean barplots summarize replicates. Heatmap log10
transform (with pseudocount) compresses the orders-of-magnitude span typical
of count matrices. `--cluster rows|cols|both` reorders the heatmap by
average-linkage hierarchical clustering (Euclidean distance on the displayed,
transformed values); pair it with `--trans log10` so abundant taxa do not
dominate the distances.

**Interpretation.** Barplots show relative structure even when the matrix
holds counts — the eye reads proportions. Compare like with like: same rank,
same top N across panels. Clustered heatmaps group taxa with similar
abundance profiles and samples with similar community structure; read blocks
of co-occurring taxa, not single rows.

**Pitfalls.** "Other" can hide the very taxon a differential test flagged —
raise `top` or plot that taxon explicitly when telling its story. Clustering
reorders but draws no dendrogram, and "Other" is clustered like any taxon;
for dendrograms, annotation bars and publication palettes hand the matrix to
the visualization plugin. Axes with fewer than two entries are left
unclustered (reported in `cluster.applied`).