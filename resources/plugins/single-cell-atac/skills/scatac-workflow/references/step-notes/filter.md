Purpose: remove low-quality cells with explicit, traceable gates. Principle:
snapatac2.filter_cells applies min_counts/max_counts on n_fragment and
min_tsse/max_tsse on TSSe in place on the backed object; each gate reports how
many cells IT removed so a surprising total traces to one threshold.
Interpretation: read the QC distributions first, then choose bounds the data
supports; a None bound skips that gate. Pitfalls: filter_cells uses a truthy
check on min_tsse, so 0 or None disables the TSSe gate — and TSSe must have been
computed (gene_anno) or the gate is skipped with a note; keeping under half the
cells (kept_fraction < 0.5) usually means a threshold is too strict, stop and
re-read the distributions. Doublet removal is NOT here — scrublet needs the
selected feature matrix and runs in the features step. Filtering is in place on
the backed file, so the before snapshot is captured as plain numbers first.