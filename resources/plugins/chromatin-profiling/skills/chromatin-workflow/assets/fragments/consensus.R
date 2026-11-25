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
if (length(all_peaks) == 0) {
  stop("no peaks remain after filtering (all input peaks may have fallen inside blacklisted regions)", call. = FALSE)
}
disjoint <- GenomicRanges::disjoin(all_peaks)
per_sample <- matrix(
  vapply(peak_list, function(peaks) {
    as.integer(IRanges::overlapsAny(disjoint, peaks))
  }, integer(length(disjoint))),
  nrow = length(disjoint),
  ncol = length(peak_list),
  dimnames = list(NULL, names(peak_list))
)
occupancy <- rowSums(per_sample)
keep <- occupancy >= effective
consensus_peaks <- GenomicRanges::reduce(disjoint[keep])
if (length(consensus_peaks) == 0) {
  stop(sprintf("no consensus peaks present in >= %d of %d samples", effective, length(peak_list)), call. = FALSE)
}

# Recompute occupancy on the merged consensus so rows align with the final peaks.
consensus_ids <- paste0(
  GenomicRanges::seqnames(consensus_peaks), ":",
  GenomicRanges::start(consensus_peaks), "-",
  GenomicRanges::end(consensus_peaks)
)
occupancy_matrix <- matrix(
  vapply(peak_list, function(peaks) {
    as.integer(IRanges::overlapsAny(consensus_peaks, peaks))
  }, integer(length(consensus_peaks))),
  nrow = length(consensus_peaks),
  ncol = length(peak_list),
  dimnames = list(consensus_ids, names(peak_list))
)

cat(sprintf("consensus: %d peaks present in >= %d of %d samples\n",
            length(consensus_peaks), effective, length(peak_list)))
cat(sprintf("union of all peaks: %d disjoint intervals, %d kept\n",
            length(disjoint), sum(keep)))
analysis_steps <- c(analysis_steps, "consensus")
