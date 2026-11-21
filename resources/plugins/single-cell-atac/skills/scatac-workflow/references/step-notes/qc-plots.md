Purpose: make the QC distributions visible so thresholds are read, not guessed.
Principle: snapatac2 renders the fragment-size distribution and TSSe enrichment
as plotly figures; the per-cell metrics (n_fragment on a log axis, frac_dup,
frac_mito) are drawn as histograms. Interpretation: look for the empty-barcode
shoulder (low n_fragment), the doublet shoulder (high n_fragment), and the TSSe
separation between signal and noise; the before/after panels (plot_tag) show what
each gate removed. Pitfalls: figures are interactive .html (plotly), not static
images — a PNG needs kaleido plus a chrome backend (static_png); without
gene_anno there is no TSSe panel and the gap is reported, not hidden. A plot is
evidence for a threshold decision, not a substitute for reading the quantiles.