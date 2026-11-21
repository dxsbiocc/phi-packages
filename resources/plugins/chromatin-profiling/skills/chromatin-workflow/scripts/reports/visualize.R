# JSON contract for chip_visualize: sample PCA on log-normalized counts and the
# peak-overlap UpSet figure. The count matrix may come from a table or from BAMs.
pca <- if (!is.null(pca_result)) {
  scores <- as.data.frame(pca_result$scores)
  scores$sample <- rownames(scores)
  list(variance_explained = pca_result$variance_explained,
       n_peaks = pca_result$n_peaks, scores = scores)
} else NULL

emit(list(
  input = tool_input,
  output = tool_output,
  n_peaks = if (!is.null(count_matrix)) nrow(count_matrix) else 0L,
  n_samples = if (!is.null(count_matrix)) ncol(count_matrix) else 0L,
  samples = chr(colnames(count_matrix)),
  normalize_method = normalize_method,
  pca = pca,
  figures = chr(unlist(figure_paths, use.names = FALSE)),
  steps = chr(analysis_steps)
))