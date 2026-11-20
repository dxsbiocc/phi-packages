# JSON contract for sc_annotate.
notes <- character(0)
if (is.null(panel_path) && is.null(mapping_path)) {
  notes <- c(notes, "no --panel and no --mapping: nothing was annotated")
}
if (!is.null(annotation_scores)) {
  unassigned <- Filter(function(entry) is.null(entry$best), annotation_scores)
  if (length(unassigned) > 0) {
    notes <- c(notes, sprintf(
      "%d cluster(s) scored below --min-score %s: %s",
      length(unassigned), min_score,
      paste(vapply(unassigned, function(entry) entry$cluster, character(1)), collapse = ", ")))
  }
  missing_genes <- unique(unlist(lapply(annotation_scores, function(entry) {
    lapply(Filter(function(s) s$genes_missing > 0, entry$scores), function(s) s$cell_type)
  })))
  if (length(missing_genes) > 0) {
    notes <- c(notes, "some panel genes are absent from this object; check the gene symbol naming (human MT- vs mouse mt-)")
  }
}

emit(list(
  input = tool_input,
  output = tool_output,
  group_by = group_by,
  label_column = if (!is.null(annotation_assignments)) label_column else NULL,
  method = annotation_method,
  panel = if (!is.null(annotation_scores)) chr(annotation_panel) else NULL,
  scores = annotation_scores,
  assignments = annotation_assignments,
  automated = annotation_automated,
  figures = chr(figure_paths[["annotation-dotplot"]]),
  notes = chr(notes)
))