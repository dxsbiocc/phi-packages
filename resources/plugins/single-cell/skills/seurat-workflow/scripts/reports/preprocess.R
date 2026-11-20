# JSON contract for sc_preprocess.
cluster_summary <- lapply(cluster_columns, function(column) {
  values <- obj[[column, drop = TRUE]]
  counts_table <- table(values)
  list(
    column = column,
    clusters = length(counts_table),
    sizes = lapply(sorted_levels(names(counts_table)), function(name) {
      list(cluster = name, cells = as.integer(counts_table[[name]]))
    })
  )
})

emit(list(
  input = tool_input,
  output = tool_output,
  assay = loaded_assay,
  steps = chr(analysis_steps),
  parameters = list(
    normalization = normalization,
    scale_factor = scale_factor,
    selection_method = selection_method,
    n_variable_features = length(VariableFeatures(obj)),
    n_pcs = length(pca_stdev),
    dims_use = dims,
    resolution = resolution,
    seed = seed
  ),
  top_variable_features = chr(head(VariableFeatures(obj), 10)),
  variance = pca_variance,
  clusters = cluster_summary,
  embeddings = chr(embeddings),
  figures = chr(c(figure_paths[["dimplot-clusters"]], figure_paths[["pca-elbow"]])),
  notes = chr(analysis_notes)
))