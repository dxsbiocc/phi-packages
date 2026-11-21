# JSON contract for chip_consensus: the merged peak universe and how many of its
# sites each sample occupies, plus blacklist removal and the overlap UpSet figure.
per_sample <- if (!is.null(occupancy_matrix)) {
  totals <- colSums(occupancy_matrix)
  lapply(names(totals), function(name) {
    list(sample = name, n_occupied = as.integer(totals[[name]]))
  })
} else NULL

emit(list(
  input = tool_input,
  output = tool_output,
  n_samples = if (!is.null(peak_list)) length(peak_list) else 0L,
  n_consensus = if (!is.null(consensus_peaks)) length(consensus_peaks) else 0L,
  min_occupancy = min_occupancy,
  occupancy_per_sample = per_sample,
  blacklist_removed = blacklist_removed,
  figures = chr(unlist(figure_paths, use.names = FALSE)),
  steps = chr(analysis_steps)
))