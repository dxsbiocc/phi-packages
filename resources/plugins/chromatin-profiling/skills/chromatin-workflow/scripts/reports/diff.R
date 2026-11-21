# JSON contract for chip_diff: differential binding from DiffBind (DESeq2/edgeR)
# or csaw (window-based, broad marks). A csaw platform gap surfaces as `error`.
result <- if (!is.null(csaw_result)) csaw_result else diff_results
if (is.null(result)) {
  emit(list(input = tool_input, output = tool_output,
            error = "no differential result was produced", steps = chr(analysis_steps)))
} else {
  n_tested <- if (!is.null(result$n_tested)) result$n_tested else result$n_regions
  emit(list(
    input = tool_input,
    output = tool_output,
    method = result$method,
    fdr_threshold = fdr_threshold,
    lfc_threshold = lfc_threshold,
    n_significant = result$n_significant,
    n_tested = n_tested,
    n_windows = result$n_windows,
    n_regions = result$n_regions,
    table = result$table,
    note = result$note,
    error = result$error,
    figures = chr(unlist(figure_paths, use.names = FALSE)),
    steps = chr(analysis_steps)
  ))
}