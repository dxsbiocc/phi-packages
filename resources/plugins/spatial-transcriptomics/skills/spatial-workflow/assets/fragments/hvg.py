# Most genes are noise across spatial units; keeping them dilutes every distance.
# HVG selection ranks genes by how much they vary beyond the mean-variance trend of
# technical noise. seurat_v3 models raw counts directly, so it runs on the counts
# layer, not on the log-normalized matrix. These HVGs drive PCA/clustering; spatially
# variable genes (SVG step) are a separate, coordinate-aware selection.
sc.pp.highly_variable_genes(adata, layer="counts", flavor=hvg_flavor, n_top_genes=n_hvg)
analysis_steps.append("variable-features")
show_plot(sc.pl.highly_variable_genes(adata, show=False), "hvg-mean-variance")
n_hvg_found = int(adata.var["highly_variable"].sum())
print(f"{n_hvg_found} highly variable genes of {adata.n_vars}; lognorm snapshot keeps all genes")