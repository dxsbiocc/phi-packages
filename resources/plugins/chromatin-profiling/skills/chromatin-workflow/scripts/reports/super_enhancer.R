# JSON contract for chip_super_enhancer.
info_obj <- if (!is.null(se_result)) se_result else list()
payload <- list(
  input = tool_input,
  tables = chr(unlist(list(info_obj$output_table), use.names = FALSE)),
  figures = chr(unlist(figure_paths[c("super-enhancer-hockeystick")], use.names = FALSE)),
  notes = chr(analysis_notes),
  steps = chr(analysis_steps)
)
for (nm in names(info_obj)) payload[[nm]] <- info_obj[[nm]]
emit(payload)
