# Read these distributions before choosing thresholds: the violin shape tells
# you where the empty-droplet shoulder and the doublet tail sit. plot_tag is
# "before" on the first pass and "after" once filtering has run.
show_plot(
    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        jitter=0.4,
        multi_panel=True,
        show=False,
    ),
    f"qc-violin-{plot_tag}",
)
show_plot(sc.pl.highest_expr_genes(adata, n_top=20, show=False), f"qc-highest-{plot_tag}")