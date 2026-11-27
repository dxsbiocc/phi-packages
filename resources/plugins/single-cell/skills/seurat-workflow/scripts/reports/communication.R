# JSON contract for sc_communication.
payload <- list(
  input = tool_input,
  output = tool_output,
  figures = chr(unlist(figure_paths[c("communication-heatmap", "communication-dotplot")], use.names = FALSE)),
  notes = c(
    chr(analysis_notes),
    "mRNA ligand-receptor co-expression indicates signaling potential, not physical cell-cell contact",
    "pass --interactions_path with an OmniPath/CellPhoneDB TSV for exhaustive receptorome screening"
  )
)
if (!is.null(ccc_info)) {
  for (nm in names(ccc_info)) payload[[nm]] <- ccc_info[[nm]]
}
emit(payload)
