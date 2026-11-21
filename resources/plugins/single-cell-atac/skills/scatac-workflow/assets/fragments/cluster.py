# Leiden cuts the expression neighbour graph built on X_spectral. resolution is the
# granularity dial — sweep it (a list) and compare cluster counts rather than trusting
# one value. A single resolution writes obs['leiden']; several write obs['leiden_<r>'].
# UMAP is for looking, not for clustering: it preserves local neighbourhoods but
# distorts distances, so the graph and Leiden never use it.
cluster_columns = []
for value in resolution:
    column = "leiden" if len(resolution) == 1 else f"leiden_{value}"
    snap.tl.leiden(adata, resolution=value, random_state=seed, key_added=column, inplace=True)
    cluster_columns.append(column)
analysis_steps.append("clusters")

if run_umap:
    snap.tl.umap(adata, use_rep="X_spectral", random_state=seed, inplace=True)
    embeddings.append("X_umap")
    show_plot(snap.pl.umap(adata, color=cluster_columns[0], show=False), "dimplot-clusters")

_sizes = obs_col(adata, cluster_columns[0])
if _sizes is not None:
    _names, _counts = np.unique(_sizes.astype(str), return_counts=True)
    print(f"{len(_names)} clusters at resolution {resolution[0]}:")
    for _n, _c in zip(_names, _counts):
        print(f"  cluster {_n}: {_c} cells")