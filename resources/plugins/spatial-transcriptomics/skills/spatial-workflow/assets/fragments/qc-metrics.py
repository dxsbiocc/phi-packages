# Three metrics decide most spatial QC calls:
#   n_genes_by_counts  genes detected per spot/bin/cell (low = empty or off-tissue unit)
#   total_counts       UMIs per unit                  (low = same)
#   pct_counts_mt      mitochondrial fraction         (high = stressed/damaged tissue)
# Cytoplasmic RNA leaks out of a torn membrane while mitochondria stay inside,
# so a high mitochondrial fraction marks damage, not biology. Spatial tissue
# tolerates more mitochondrial signal than dissociated cells, so mt_max is looser.
# Ribosomal and hemoglobin fractions are context: they vary by region and contamination.
# mt_pattern comes from the config cell: "^MT-" for human, "^mt-" for mouse.
adata.var["mt"] = adata.var_names.str.contains(mt_pattern, regex=True)
adata.var["ribo"] = adata.var_names.str.contains(ribo_pattern, regex=True)
adata.var["hb"] = adata.var_names.str.contains(hb_pattern, regex=True)
qc_vars = [name for name in ("mt", "ribo", "hb") if adata.var[name].any()]
# percent_top defaults to 500, which raises on objects with fewer genes
sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=qc_vars,
    inplace=True,
    log1p=False,
    percent_top=tuple(p for p in (50, 100, 200, 500) if p < adata.n_vars) or None,
)
qc_scored = ["pct_counts_" + name for name in qc_vars]
obj_before = adata

# The quantiles are what a threshold should be based on, not a round number.
# pct_counts_mt only exists when mt_pattern matched genes (mouse uses ^mt-).
_qc_cols = [c for c in ("n_genes_by_counts", "total_counts", "pct_counts_mt") if c in adata.obs]
print(adata.obs[_qc_cols].describe(percentiles=[0.01, 0.5, 0.99]))