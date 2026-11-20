# Marker genes are what turn a cluster number into a cell type. FindAllMarkers
# tests every cluster against all the others. The p-values rank genes; they are
# not clean inference, because the clusters were defined on this same data.
if ("SCT" %in% names(obj@assays) && DefaultAssay(obj) == "SCT") {
  # FindAllMarkers on SCT needs the corrected counts prepared first.
  obj <- tryCatch(PrepSCTFindMarkers(obj, verbose = FALSE), error = function(e) obj)
}

markers <- FindAllMarkers(
  obj,
  group.by = group_by,
  test.use = test_use,
  only.pos = only_positive,   # a marker is enriched in its cluster
  min.pct = min_pct,          # must be detected in at least this fraction
  logfc.threshold = logfc_threshold,
  max.cells.per.ident = if (is.null(max_cells_per_ident)) Inf else max_cells_per_ident,
  verbose = FALSE
)

if (nrow(markers) == 0) {
  message("no marker passed the thresholds; lower logfc_threshold or min_pct, or check that clusters differ at all")
} else {
  markers$cluster <- as.character(markers$cluster)
  # Keep the full table on the object so annotate and downstream tools reuse it.
  obj@misc[["markers"]] <- markers

  # One machine-readable table per cluster, sorted by significance then effect.
  marker_tables <- character(0)
  for (group in sort(unique(markers$cluster))) {
    rows <- markers[markers$cluster == group, , drop = FALSE]
    rows <- rows[order(rows$p_val_adj, -rows$avg_log2FC), , drop = FALSE]
    if (!is.null(tables_dir) && nzchar(tables_dir)) {
      dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
      path <- file.path(tables_dir, paste0("markers-", gsub("[^A-Za-z0-9._-]", "_", group), ".csv"))
      write.csv(rows, path, row.names = FALSE)
      marker_tables <- c(marker_tables, path)
    }
    print(head(rows[, c("gene", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")], top_n))
  }

  # avg_log2FC is a fold change from this test only, not an effect size you can
  # compare across datasets. pct.1/pct.2 (detection fractions) are the specificity.
  if (plot_heatmap) {
    top_genes <- unique(unlist(lapply(sort(unique(markers$cluster)), function(group) {
      rows <- markers[markers$cluster == group, , drop = FALSE]
      rows <- rows[order(rows$p_val_adj, -rows$avg_log2FC), , drop = FALSE]
      head(rows$gene, top_n)
    })))
    show_plot(DoHeatmap(obj, features = head(top_genes, 100), group.by = group_by, size = 3) + NoLegend(),
              paste0("markers-heatmap-", gsub("[^A-Za-z0-9._-]", "_", group_by)),
              width = 10, height = 8)
  }
}