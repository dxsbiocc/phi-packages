# Neighbours are built in PC space, never on the 2D embedding: UMAP/t-SNE keep
# local neighbourhoods but distort distances and densities, so clustering on them
# invents structure. Leiden then cuts the neighbour graph; resolution is the
# granularity dial — sweep it and compare, do not trust one value.
use_rep = "X_pca_integrated" if "X_pca_integrated" in adata.obsm else "X_pca"
if use_rep == "X_pca_integrated":
    analysis_notes.append("neighbours use the integrated embedding X_pca_integrated")
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=dims, use=f"obsm.{use_rep}", rng=seed)
analysis_steps.append("neighbors")

cluster_columns = []
for value in resolution:
    column = "leiden" if len(resolution) == 1 else f"leiden_{value}"
    sc.tl.leiden(adata, resolution=value, rng=seed, out=f"obs.{column}", flavor="igraph")
    cluster_columns.append(column)
analysis_steps.append("clusters")

embeddings = []
if run_umap:
    sc.tl.umap(adata, rng=seed)
    embeddings.append("umap")
if run_tsne:
    sc.tl.tsne(adata, rng=seed, use=f"obsm.{use_rep}")
    embeddings.append("tsne")

if embeddings:
    plot = (
        sc.pl.umap(adata, color=cluster_columns[0], show=False)
        if embeddings[0] == "umap"
        else sc.pl.tsne(adata, color=cluster_columns[0], show=False)
    )
    show_plot(plot, "dimplot-clusters")
print(adata.obs[cluster_columns[0]].value_counts().sort_index())
if dims > n_pcs:
    analysis_notes.append(f"dims_use ({dims}) exceeds computed PCs ({n_pcs})")