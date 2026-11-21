# A consensus peakset: the union of reproducible binding sites across samples.
# Disjoin every sample's peaks into non-overlapping intervals, keep intervals
# present in at least min_occupancy samples, then merge adjacent kept intervals.
# This is the region universe the count matrix and differential analysis use, so
# every sample is measured over the same coordinates.
if (is.null(peak_list) || length(peak_list) == 0) {
  stop("load peaks before building a consensus (load-peaks)", call. = FALSE)
}
effective <- min(min_occupancy, length(peak_list))
if (effective < 1) effective <- 1L

all_peaks <- unlist(peak_list, use.names = FALSE)
disjoint <- GenomicRanges::disjoin(all_peaks)
per_sample <- vapply(peak_list, function(peaks) {
  as.integer(IRanges::overlapsAny(disjoint, peaks))
}, integer(length(disjoint)))
occupancy <- rowSums(per_sample)
keep <- occupancy >= effective
consensus_peaks <- GenomicRanges::reduce(disjoint[keep])

# Recompute occupancy on the merged consensus so rows align with the final peaks.
occupancy_matrix <- vapply(peak_list, function(peaks) {
  as.integer(IRanges::overlapsAny(consensus_peaks, peaks))
}, integer(length(consensus_peaks)))
rownames(occupancy_matrix) <- paste0(seqnames(consensus_peaks), ":",
                                     start(consensus_peaks), "-", end(consensus_peaks))

cat(sprintf("consensus: %d peaks present in >= %d of %d samples\n",
            length(consensus_peaks), effective, length(peak_list)))
cat(sprintf("union of all peaks: %d disjoint intervals, %d kept\n",
            length(disjoint), sum(keep)))
analysis_steps <- c(analysis_steps, "consensus")