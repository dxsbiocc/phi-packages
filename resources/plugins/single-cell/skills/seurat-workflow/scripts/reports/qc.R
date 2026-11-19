# JSON contract for sc_qc. Reads the bookkeeping the fragments left behind.
before <- list(cells = ncol(obj_before), genes = nrow(obj_before), per_sample = NULL)
if (sample_column %in% colnames(obj_before[[]])) {
  counts_table <- table(obj_before[[sample_column, drop = TRUE]])
  before$per_sample <- lapply(names(counts_table), function(name) {
    list(sample = name, cells = as.integer(counts_table[[name]]))
  })
}

metrics <- obj_before[[]][, c("nFeature_RNA", "nCount_RNA", "percent.mt"), drop = FALSE]
distribution <- lapply(names(metrics), function(column) {
  values <- metrics[[column]]
  list(
    metric = column,
    min = round(min(values, na.rm = TRUE), 3),
    q01 = round(stats::quantile(values, quantiles[[1]], na.rm = TRUE), 3),
    median = round(stats::median(values, na.rm = TRUE), 3),
    q99 = round(stats::quantile(values, quantiles[[length(quantiles)]], na.rm = TRUE), 3),
    max = round(max(values, na.rm = TRUE), 3)
  )
})

after <- list(
  cells = ncol(obj),
  genes = nrow(obj),
  removed_cells = before$cells - ncol(obj),
  removed_genes = before$genes - nrow(obj),
  kept_fraction = round(ncol(obj) / max(1, before$cells), 4)
)

warnings <- character(0)
if (ncol(obj) == 0) warnings <- c(warnings, "all cells were removed; relax the thresholds")
if (after$kept_fraction < 0.5 && ncol(obj) > 0) {
  warnings <- c(warnings, sprintf(
    "only %.1f%% of cells survived; check whether a threshold is too strict", after$kept_fraction * 100))
}
if (!is.null(doublet_result$error) && detect_doublets) {
  warnings <- c(warnings, paste0("doublet detection unavailable: ", doublet_result$error))
}

emit(list(
  input = tool_input,
  output = tool_output,
  before = before,
  after = after,
  distribution = distribution,
  filters = qc_filters,
  features_scored = chr(qc_scored),
  doublets = doublet_result,
  figures = chr(c(figure_paths[["qc-violin-before"]], figure_paths[["qc-violin-after"]])),
  warnings = chr(warnings)
))