# PCA rotates the variable-gene space into orthogonal axes ordered by variance.
# The mask keeps all genes in the object (markers need them) while the axes are
# computed from variable genes only. Read the elbow: axes past the bend are
# mostly noise, but keep a few extra so rare populations survive.
sc.tl.pca(adata, n_comps=n_pcs, mask="var.highly_variable", svd_solver="arpack", rng=seed)
analysis_steps.append("pca")
variance_ratio = adata.uns["pca"]["variance_ratio"]
share = [round(float(v) * 100, 3) for v in variance_ratio]
pca_variance = {
    "pcs": len(share),
    "variance_explained": share,
    "first_ten_percent": round(sum(share[:10]), 3),
}
show_plot(sc.pl.pca_variance_ratio(adata, log=False, show=False), "pca-elbow")
print("variance explained by first PCs:", np.round(variance_ratio[:10], 3))