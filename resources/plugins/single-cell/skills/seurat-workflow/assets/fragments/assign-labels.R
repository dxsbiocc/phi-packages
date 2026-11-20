# Write the labels back. Edit this by hand if your reading of the dot plot
# disagrees with the automatic call -- that judgement is the analysis.
# A cluster with no score stays unlabelled rather than inheriting a wrong name.
groups <- sort(unique(as.character(obj[[group_by, drop = TRUE]])))

if (!is.null(mapping_path) && nzchar(mapping_path)) {
  if (!file.exists(mapping_path)) stop("mapping not found: ", mapping_path, call. = FALSE)
  mapping <- unlist(jsonlite::fromJSON(mapping_path, simplifyVector = TRUE))
  names(mapping) <- gsub("-", "_", names(mapping))
  missing_clusters <- setdiff(groups, names(mapping))
  if (length(missing_clusters) > 0) {
    stop("mapping is missing clusters: ", paste(missing_clusters, collapse = ", "),
         "; label every cluster or omit mapping_path", call. = FALSE)
  }
  labels <- unname(mapping[as.character(obj[[group_by, drop = TRUE]])])
  annotation_assignments <- lapply(groups, function(cluster) {
    list(cluster = cluster, label = unname(mapping[[cluster]]))
  })
} else {
  best_labels <- vapply(annotation_scores, function(e) {
    if (is.null(e$best)) "unassigned" else e$best$cell_type
  }, character(1))
  names(best_labels) <- vapply(annotation_scores, function(e) e$cluster, character(1))
  if (all(best_labels == "unassigned")) {
    stop("No cluster could be scored against the panel; fix the identifiers before labelling.",
         call. = FALSE)
  }
  labels <- unname(best_labels[as.character(obj[[group_by, drop = TRUE]])])
  annotation_assignments <- lapply(annotation_scores, function(e) {
    list(
      cluster = e$cluster,
      label = if (is.null(e$best)) "unassigned" else e$best$cell_type,
      score = if (is.null(e$best)) NA_real_ else e$best$score
    )
  })
}

obj[[label_column]] <- labels
table(obj[[label_column, drop = TRUE]])

if ("umap" %in% names(obj@reductions)) {
  show_plot(DimPlot(obj, reduction = "umap", group.by = label_column, label = TRUE),
            "dimplot-celltypes", width = 7, height = 6)
}