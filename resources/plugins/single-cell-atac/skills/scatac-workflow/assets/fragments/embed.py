# Spectral embedding is the scATAC analogue of LSI/PCA: a diffusion-map embedding
# of the cell graph built on the selected tile/peak features (cosine distance,
# weighted by feature standard deviation). It denoises the sparse binary accessibility
# matrix into a low-dimensional space where neighbours are meaningful. Unlike PCA it
# has no variance-ratio curve — choose n_comps from the downstream clustering, not a
# scree plot. The kNN graph is built on X_spectral (never on UMAP, which distorts
# distances); Leiden cuts this graph in the cluster step.
snap.tl.spectral(
    adata, n_comps=n_comps, features="selected", random_state=seed,
    num_threads=max(1, n_jobs), inplace=True,
)
analysis_steps.append("spectral")
snap.pp.knn(adata, n_neighbors=n_neighbors, use_rep="X_spectral", random_state=seed, inplace=True)
analysis_steps.append("knn")
embeddings.append("X_spectral")
print(f"spectral embedding X_spectral ({n_comps} comps); kNN graph k={n_neighbors}")