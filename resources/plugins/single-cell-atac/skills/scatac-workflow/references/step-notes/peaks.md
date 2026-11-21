Purpose: turn uniform tiles into the open-chromatin regions that differ between
cell states. Principle: snapatac2.tl.macs3 calls peaks per cluster on that
cluster's fragments (MACS3 models the fragment pileup against a local background);
merge_peaks extends each summit by peak_half_width and merges overlaps into one
non-overlapping set; make_peak_matrix builds a SEPARATE peak-count object
(peak_adata) with the cluster labels copied across. Interpretation: peaks per
cluster and the merged peak count tell you how much regulatory structure the
data supports; the merged BED is the feature set of record. Pitfalls: MACS3 needs
fragments — a peak-matrix input without fragments reports a gap; the backed tile
matrix cannot be resized in place, so the peak matrix is a new file (inplace=False);
call_broad_peaks is for histone marks, ATAC open chromatin uses narrow peaks;
peak_qvalue and peak_half_width are biological decisions about sensitivity and
region width. Peaks drive markers; tiles stay the embedding.