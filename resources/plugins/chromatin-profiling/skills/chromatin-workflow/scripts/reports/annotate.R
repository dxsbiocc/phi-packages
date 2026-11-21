# JSON contract for chip_annotate: genomic feature distribution, TSS region and
# the annotation source. Handles a single csAnno or a per-sample list of them.
anno_df <- if (inherits(annotation_result, "csAnno")) {
  as.data.frame(annotation_result)
} else if (is.list(annotation_result)) {
  do.call(rbind, lapply(annotation_result, as.data.frame))
} else NULL

feature_distribution <- if (!is.null(anno_df) && "annotation" %in% colnames(anno_df)) {
  counts <- table(anno_df$annotation)
  lapply(names(counts), function(feature) {
    list(feature = feature, n_peaks = as.integer(counts[[feature]]))
  })
} else NULL

source <- if (!is.null(txdb_path) && nzchar(txdb_path)) {
  txdb_path
} else if (!is.null(gtf_path) && nzchar(gtf_path)) {
  gtf_path
} else {
  NULL
}

emit(list(
  input = tool_input,
  output = tool_output,
  n_peaks = if (!is.null(anno_df)) nrow(anno_df) else 0L,
  annotation_source = source,
  tss_region = c(-tss_upstream, tss_downstream),
  feature_distribution = feature_distribution,
  gene_symbols = !is.null(org_db) && nzchar(org_db) && have(org_db),
  figures = chr(unlist(figure_paths, use.names = FALSE)),
  steps = chr(analysis_steps)
))