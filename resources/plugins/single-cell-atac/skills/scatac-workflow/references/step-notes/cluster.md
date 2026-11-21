Purpose: partition cells into clusters and lay them out for visualization.
Principle: snapatac2.tl.leiden optimizes modularity on the kNN graph at a given
resolution; UMAP embeds X_spectral into 2D for display. Interpretation:
resolution controls granularity — sweep several resolutions (a list) to compare
stable vs over-split partitions; read cluster sizes, a cluster of a handful of
cells may be doublets or noise. Pitfalls: clusters are defined from the same
data used to test markers (double dipping), so marker p-values are descriptive;
the cluster column (leiden, or leiden_<resolution> for a sweep) is what peak
calling, markers and annotation group by — keep group_by consistent; UMAP is a
visualization, distances on it are not quantitative. seed makes clustering and
embedding repeat exactly (snapatac2 default 0).