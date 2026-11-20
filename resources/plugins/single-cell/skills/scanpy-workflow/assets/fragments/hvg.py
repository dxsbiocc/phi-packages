# Most genes are noise across cells; keeping them dilutes every distance. HVG
# selection ranks genes by how much they vary beyond the mean-variance trend of
# technical noise. seurat_v3 models raw counts directly, so it runs on the counts
# layer, not on the log-normalized matrix.
hvg_kwargs = {"use": "layers.counts", "flavor": hvg_flavor, "n_top_genes": n_hvg}
if batch_key and batch_key in adata.obs.columns:
    # Per-batch ranking stops one dominant sample from choosing every gene.
    hvg_kwargs["batch"] = f"obs.{batch_key}"
elif batch_key:
    analysis_notes.append(f"batch_key '{batch_key}' not in obs; HVG ran without batches")
sc.pp.highly_variable_genes(adata, **hvg_kwargs)
analysis_steps.append("variable-features")
show_plot(sc.pl.highly_variable_genes(adata, show=False), "hvg-mean-variance")
n_hvg_found = int(adata.var["highly_variable"].sum())
print(f"{n_hvg_found} highly variable genes of {adata.n_vars}; lognorm snapshot keeps all genes")