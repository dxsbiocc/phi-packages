Purpose: differential binding with DiffBind + DESeq2 — the standard route for TF
and narrow-mark ChIP, and for CUT&Tag. Principle: DiffBind builds a consensus count
matrix, normalizes for library size, and runs DESeq2's negative-binomial test per
peak; a contrast compares two conditions, and with no explicit contrast every
condition pair with enough replicates is tested. Interpretation: dba.contrast
(contrast=c(col,treated,control) or categories=DBA_CONDITION, minMembers) sets the
comparison; dba.analyze(DBA_DESEQ2) fits; dba.report(th=fdr, fold=lfc,
bUsePval=FALSE) returns FDR-sorted significant peaks. n_significant of n_tested is
the headline; the MA plot shows fold-change vs abundance. Pitfalls: the sample is
the experimental unit — fewer than min_members (default 2) replicates per
condition means no condition-level claim. DESeq2 is more conservative than edgeR.
Zero significant peaks with a clean run is a null or under-powered result, not "no
regulation". Runs on raw counts with internal normalization, never the normalized
exploration matrix.