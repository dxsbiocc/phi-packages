# JSON contract for sc_pseudobulk.
emit(list(
  input = tool_input,
  outputs = chr(pb_info$written),
  assay = loaded_assay,
  units = ncol(pb_info$matrix),
  genes = nrow(pb_info$matrix),
  design = pb_info$design,
  sample_column = sample_column,
  group_column = if (pb_info$group_provided) group_column else NULL,
  wrapper_hint = if (pb_info$group_provided) {
    list(
      counts = counts_output,
      samplesheet = samples_output,
      sample_id_col = "sample",
      gene_id_col = "gene_id",
      contrast_variable = pb_info$condition_column
    )
  } else {
    NULL
  },
  parameters = list(
    min_cells = min_cells,
    min_count = min_count,
    min_cells_expressing = min_cells_expressing
  ),
  dropped_genes = pb_info$dropped,
  notes = chr(pb_info$notes)
))