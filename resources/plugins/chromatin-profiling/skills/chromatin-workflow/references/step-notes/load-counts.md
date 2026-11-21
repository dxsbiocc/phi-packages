Purpose: read an existing peak-by-sample count matrix — the alternative entry
point when counting already happened upstream (featureCounts, or a DiffBind
dba.count export). Principle: rows are peaks, columns are samples, the first
column holds peak ids and the header holds sample names; the matrix is coerced to
numeric and rejected if any entry is non-numeric. Interpretation: it reports
n_peaks, n_samples and per-sample library sizes, then feeds normalize/pca/diff
directly without re-counting. Pitfalls: a count-table input carries no peaks or
BAM, so FRiP, consensus and coverage profiles are unavailable from it — say which
input you analysed. The separator is comma for .csv, tab otherwise; a wrong
separator silently mis-parses the matrix. These counts are raw and may feed a
differential engine; the normalized copy is for exploration only.