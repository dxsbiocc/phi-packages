Purpose: assess whether biological replicates of the same condition call the same
peaks. Principle: the formal test is IDR (irreproducible discovery rate), which
needs the `idr` pip package this environment does not ship — a reported gap. This
step gives a transparent proxy instead: the fraction of each sample's peaks that
overlap another sample (GenomicRanges::overlapsAny), plus a Spearman count
correlation when a count matrix exists. Interpretation: same-condition replicate
pairs (from the sheet's Condition column) are the overlap that matters most; high
overlap means reproducible binding, near-zero means biological or technical
divergence. Pitfalls: overlap fraction is NOT IDR — never present it as a
calibrated reproducibility test; it needs at least two peak sets or it skips. The
proxy is descriptive evidence to flag failed replicates, not a significance
statement.