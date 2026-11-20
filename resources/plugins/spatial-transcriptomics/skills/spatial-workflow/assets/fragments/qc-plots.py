# Read these distributions before choosing thresholds. In spatial data the low
# shoulder is empty / off-tissue spots or bins, not empty droplets; the high tail
# can be a merged bin or a region of dense tissue. plot_tag is "before" on the
# first pass and "after" once filtering has run.
qc_columns = ["n_genes_by_counts", "total_counts"] + [c for c in qc_scored if c in adata.obs]
show_plot(
    sc.pl.violin(adata, qc_columns, jitter=0.4, multi_panel=True, show=False),
    f"qc-violin-{plot_tag}",
)
show_plot(sc.pl.highest_expr_genes(adata, n_top=20, show=False), f"qc-highest-{plot_tag}")

# The same metrics mapped back onto tissue: a QC gradient that follows anatomy
# (e.g. low counts at the tissue edge) is biology, not a threshold to filter on.
if coord_key in adata.obsm:
    show_plot(
        sc.pl.embedding(adata, basis=coord_key, color=qc_columns, show=False),
        f"qc-spatial-{plot_tag}",
    )