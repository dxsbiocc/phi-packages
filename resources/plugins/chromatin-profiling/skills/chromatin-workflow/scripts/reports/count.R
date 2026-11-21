# JSON contract for chip_count: the consensus peak x sample affinity matrix, its
# library sizes and the normalization applied for exploration.
library_sizes <- if (!is.null(count_matrix)) {
  totals <- colSums(count_matrix)
  lapply(names(totals), function(name) {
    list(sample = name, total_reads = as.integer(totals[[name]]))
  })
} else NULL

emit(list(
  input = tool_input,
  output = tool_output,
  n_peaks = if (!is.null(count_matrix)) nrow(count_matrix) else 0L,
  n_samples = if (!is.null(count_matrix)) ncol(count_matrix) else 0L,
  samples = chr(colnames(count_matrix)),
  library_sizes = library_sizes,
  normalize_method = normalize_method,
  figures = chr(unlist(figure_paths, use.names = FALSE)),
  steps = chr(analysis_steps)
))