Purpose: library-size normalization so samples are comparable for exploration.
Principle: CPM (counts per million, the libsize default) removes sequencing-depth
differences; RPKM additionally divides by peak width, which matters when comparing
peaks of very different lengths; none keeps raw counts. When a DBA object exists,
DiffBind::dba.normalize(method=DBA_NORM_LIB) normalizes it for dba.analyze.
Interpretation: the normalized matrix feeds PCA and occupancy visualization. A
sample with zero total counts cannot normalize meaningfully (a warning is
emitted). Pitfalls: this matrix is for exploration, NOT the differential test —
DESeq2/edgeR/csaw do their own internal normalization (median-of-ratios/TMM) on raw
counts, and feeding them external CPM breaks the mean-variance relationship they
assume. RPKM needs consensus peak widths, so it requires a consensus (not a bare
count table).