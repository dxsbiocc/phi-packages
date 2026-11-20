# The integration verdict is visual: before, samples separate; after, they mix
# within each cluster. If a whole biological population disappears after
# correction, the method over-corrected — go back, do not proceed.
if batch_key not in adata.obs.columns or "X_pca_integrated" not in adata.obsm:
    print("Nothing to compare: no batch column or no corrected representation.")
else:
    sc.pp.neighbors(adata, n_pcs=dims, use="obsm.X_pca", rng=seed, out="obsp.uncorrected")
    sc.tl.umap(adata, use="obsp.uncorrected", out="obsm.X_umap_uncorrected", rng=seed)
    show_plot(
        sc.pl.embedding(adata, basis="umap_uncorrected", color=batch_key, title="before integration", show=False),
        "batch-before",
    )
    sc.pp.neighbors(adata, n_pcs=dims, use="obsm.X_pca_integrated", rng=seed)
    sc.tl.umap(adata, rng=seed)
    show_plot(sc.pl.umap(adata, color=batch_key, title="after integration", show=False), "batch-after")