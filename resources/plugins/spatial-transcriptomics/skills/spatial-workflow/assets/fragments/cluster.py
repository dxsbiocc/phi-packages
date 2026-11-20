# These neighbours are built in EXPRESSION PC space, never on the 2D embedding and
# never on tissue coordinates: UMAP keeps local neighbourhoods but distorts distances,
# so clustering on it invents structure. The spatial adjacency graph used by SVG /
# neighbourhood / communication is built separately in the spatial-neighbors step.
# Leiden cuts the expression neighbour graph; resolution is the granularity dial —
# sweep it and compare, do not trust one value.
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=dims, use_rep="X_pca", random_state=seed)
analysis_steps.append("neighbors")

cluster_columns = []
for value in resolution:
    column = "leiden" if len(resolution) == 1 else f"leiden_{value}"
    sc.tl.leiden(adata, resolution=value, random_state=seed, key_added=column, flavor="igraph")
    cluster_columns.append(column)
analysis_steps.append("clusters")

embeddings = []
if run_umap:
    sc.tl.umap(adata, random_state=seed)
    embeddings.append("umap")
    show_plot(sc.pl.umap(adata, color=cluster_columns[0], show=False), "dimplot-clusters")
print(adata.obs[cluster_columns[0]].value_counts().sort_index())
if dims > n_pcs:
    analysis_notes.append(f"dims ({dims}) exceeds computed PCs ({n_pcs})")