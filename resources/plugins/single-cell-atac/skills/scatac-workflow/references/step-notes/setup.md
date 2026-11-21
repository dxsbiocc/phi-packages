Purpose: set up the shared namespace every later step reads — imports, helpers,
bookkeeping globals and the figure writer. Principle: snapatac2 is a Rust engine
over a file-backed AnnData; the helpers (obs_has/obs_col/var_col/sorted_levels)
exist because a backed object's obs/var are lazy Rust collections, not pandas
frames, so `.columns`/`.keys()`/iteration can panic while item access and `in`
work. Interpretation: nothing biological happens here; it defines the contract
the fragments and the JSON reports share. Pitfalls: a backed object must be
closed before exit or the HDF5 corrupts; snapatac2 figures are plotly Figures
(written as interactive .html), not matplotlib Axes — only the gene-activity
downstream uses scanpy/matplotlib. The genome's chrom_sizes are bundled for
hg38/mm10/hg19/mm39; a custom genome needs a chrom_sizes file.