Purpose: turn the neighbour graph into discrete hypotheses about cell states.
Principle: neighbours are built in PC space (corrected embedding when present),
then Leiden cuts the graph into communities; resolution is the granularity dial.
Interpretation: sweep resolutions and compare (clustree, ARI, silhouette):
stable clusters across resolutions are trustworthy, and cluster size matters —
12 cells and 12,000 deserve different scrutiny. Clusters are not cell types; a
cluster earns a name only from markers. Pitfalls: never cluster on the 2D
embedding — UMAP/t-SNE preserve local neighbourhoods but distort distances and
densities; subcluster with fresh HVG/PCA, not by re-cutting the old graph.
