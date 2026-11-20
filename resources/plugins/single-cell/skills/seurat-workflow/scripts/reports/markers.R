# JSON contract for sc_markers.
if (is.null(markers) || nrow(markers) == 0) {
  emit(list(
    input = tool_input,
    output = tool_output,
    group_by = group_by,
    test = test_use,
    total_markers = 0L,
    groups = list(),
    notes = chr(c(
      "no marker passed the thresholds; lower --logfc-threshold or --min-pct, ",
      "or check that clusters differ at all"
    ))
  ))
  quit(status = 0, save = "no")
}

groups <- sorted_levels(markers$cluster)
per_group <- lapply(groups, function(group) {
  rows <- markers[markers$cluster == group, , drop = FALSE]
  rows <- rows[order(rows$p_val_adj, -rows$avg_log2FC), , drop = FALSE]
  top <- head(rows, top_n)
  list(
    cluster = group,
    markers = nrow(rows),
    top_genes = chr(head(top$gene, top_n)),
    top_log2fc = lapply(seq_len(nrow(top)), function(index) {
      list(gene = top$gene[[index]], avg_log2FC = round(top$avg_log2FC[[index]], 4))
    })
  )
})

notes <- character(0)
if (only_positive) notes <- c(notes, "only positive markers are reported (--only-positive)")
if (!"presto" %in% rownames(utils::installed.packages())) {
  notes <- c(notes, "presto is not installed; the Wilcoxon test falls back to Seurat's own implementation and is slower")
}

emit(list(
  input = tool_input,
  output = tool_output,
  group_by = group_by,
  test = test_use,
  total_markers = nrow(markers),
  groups = per_group,
  tables = chr(marker_tables),
  figures = chr(figure_paths[[paste0("markers-heatmap-", gsub("[^A-Za-z0-9._-]", "_", group_by))]]),
  notes = chr(notes)
))