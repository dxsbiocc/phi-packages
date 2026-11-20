# Annotation maps clusters to cell types. Two routes:
#   1. score each cluster against a panel of known markers (shown here)
#   2. assign labels from literature knowledge, then verify with the dot plot
# The panel is a named list (panel_path JSON, or the example below). Genes absent
# from your data are reported, not silently ignored.
if (!is.null(panel_path) && nzchar(panel_path)) {
  if (!file.exists(panel_path)) stop("panel not found: ", panel_path, call. = FALSE)
  parsed <- jsonlite::fromJSON(panel_path, simplifyVector = FALSE)
  if (!is.list(parsed) || is.null(names(parsed))) {
    stop("panel JSON must be an object mapping cell type -> array of gene symbols", call. = FALSE)
  }
  panel <- lapply(parsed, function(genes) as.character(unlist(genes)))
} else {
  panel <- list(
    T_cell      = c("CD3D", "CD3E", "TRAC", "IL7R", "LTB"),
    B_cell      = c("MS4A1", "CD79A", "CD79B", "CD19"),
    Plasma      = c("MZB1", "XBP1", "JCHAIN", "SDC1"),
    NK          = c("NKG7", "KLRD1", "GNLY", "PRF1"),
    Myeloid     = c("LYZ", "S100A8", "S100A9", "CTSS"),
    Dendritic   = c("FCER1A", "CLEC10A", "CD1C", "ITGAX"),
    Fibroblast  = c("COL3A1", "PDGFRA", "DCN", "LAMA2"),
    Epithelial  = c("EPCAM", "KRT8", "KRT18", "CDH1"),
    Endothelial = c("PECAM1", "VWF", "CDH5", "CLDN5")
  )
}
annotation_panel <- names(panel)

# Which panel genes this dataset actually contains, and which are missing.
present <- unique(unlist(panel))
missing <- setdiff(present, rownames(obj))
cat("panel genes missing from this object:", length(missing), "\n")
if (length(missing) > 0) print(missing)

# An entirely absent panel means the gene identifiers do not match this object
# (wrong species, Ensembl IDs instead of symbols). Say so and stop, rather than
# emitting labels that came from no evidence.
usable <- unique(unlist(lapply(panel, function(g) intersect(g, rownames(obj)))))
if (length(usable) == 0) {
  stop(paste(
    "None of the panel genes occur in this object.",
    "Check the identifier type: rownames(obj) holds",
    paste(head(rownames(obj), 3), collapse = ", "), "...",
    "Human symbols use MT-, mouse mt-. Ensembl IDs need a different panel."
  ), call. = FALSE)
}

# Score = mean expression of the panel genes in the cluster minus the mean of
# all other genes, so panels of different sizes stay comparable. A narrow win
# between two similar panels is not an identification: read the margin.
expression <- LayerData(obj, layer = "data")
all_genes <- rownames(expression)
groups <- sort(unique(as.character(obj[[group_by, drop = TRUE]])))

annotation_scores <- lapply(groups, function(group) {
  cells <- colnames(obj)[obj[[group_by, drop = TRUE]] == group]
  group_mean <- Matrix::rowMeans(expression[, cells, drop = FALSE])
  scores <- lapply(names(panel), function(cell_type) {
    genes <- intersect(panel[[cell_type]], all_genes)
    if (length(genes) == 0) {
      return(list(cell_type = cell_type, score = NA_real_, genes_found = 0L,
                  genes_missing = length(panel[[cell_type]])))
    }
    rest <- setdiff(all_genes, genes)
    panel_mean <- mean(group_mean[genes])
    rest_mean <- if (length(rest) > 0) mean(group_mean[rest]) else 0
    list(
      cell_type = cell_type,
      score = round(panel_mean - rest_mean, 4),
      panel_mean = round(panel_mean, 4),
      genes_found = length(genes),
      genes_missing = length(panel[[cell_type]]) - length(genes)
    )
  })
  scored <- Filter(function(entry) !is.na(entry$score), scores)
  scored <- scored[order(-vapply(scored, function(entry) entry$score, numeric(1)))]
  best <- if (length(scored) > 0 && scored[[1]]$score >= min_score) scored[[1]] else NULL
  list(cluster = group, cells = length(cells), best = best, scores = scores)
})

best_table <- data.frame(
  cluster = groups,
  label = vapply(annotation_scores, function(e) if (is.null(e$best)) NA_character_ else e$best$cell_type, character(1)),
  score = vapply(annotation_scores, function(e) if (is.null(e$best)) NA_real_ else e$best$score, numeric(1)),
  runner_up = vapply(annotation_scores, function(e) {
    scored <- Filter(function(s) !is.na(s$score), e$scores)
    if (length(scored) < 2) NA_character_ else scored[[2]]$cell_type
  }, character(1)),
  margin = vapply(annotation_scores, function(e) {
    scored <- Filter(function(s) !is.na(s$score), e$scores)
    if (length(scored) < 2) NA_real_ else round(scored[[1]]$score - scored[[2]]$score, 4)
  }, numeric(1))
)
print(best_table)

# Reference-based annotation is an alternative when a matching reference exists;
# it still needs marker verification, and it forces a label even for populations
# the reference never saw.
annotation_automated <- NULL
if (identical(annotation_method, "singler")) {
  if (!have("SingleR") || !have("celldex")) {
    annotation_automated <- list(
      method = "SingleR",
      error = paste(
        "SingleR and celldex are not installed in the plugin environment;",
        "install them into a project environment, or use --method panel with a marker panel"
      )
    )
    message(annotation_automated$error)
  } else {
    annotation_automated <- tryCatch({
      reference <- celldex::HumanPrimaryCellAtlasData()
      prediction <- SingleR::SingleR(test = expression, ref = reference, labels = reference$label.main)
      labels <- prediction$labels
      names(labels) <- colnames(obj)
      obj[["singler_label"]] <- labels
      counts <- table(labels)
      list(
        method = "SingleR",
        reference = "celldex::HumanPrimaryCellAtlasData",
        per_cell = length(labels),
        distinct = length(counts),
        top = lapply(sort(names(counts)), function(name) list(label = name, cells = as.integer(counts[[name]])))
      )
    }, error = function(error) list(method = "SingleR", error = conditionMessage(error)))
  }
}

# Only plot genes that exist: DotPlot errors on any absent feature.
plot_genes <- intersect(head(unique(unlist(panel)), 60), rownames(obj))
if (length(plot_genes) > 0) {
  show_plot(DotPlot(obj, features = plot_genes, group.by = group_by) +
              ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)),
            "annotation-dotplot", width = 12, height = 6)
}