Purpose: group spatial units into domains in expression space and embed them.
Principle: the neighbor graph is built in expression PC space (never on the 2D
UMAP and never on tissue coordinates — UMAP distorts distances, so clustering on
it invents structure); Leiden cuts that graph, and resolution is the granularity
dial. Interpretation: sweep resolution and compare; a cluster is a spatial
domain only when its units sit together in tissue (the spatial-plots step maps
them back). A scattered cluster is a cell state spread across regions.
Pitfalls: trusting one resolution; the expression neighbor graph here is
DIFFERENT from the spatial adjacency graph built in spatial-neighbors — SVG and
neighbourhood statistics use the latter, not this one. Cluster IDs are
hypotheses until markers say what they are.