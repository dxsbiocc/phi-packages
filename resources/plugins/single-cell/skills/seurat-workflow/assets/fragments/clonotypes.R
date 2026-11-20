# Cell Ranger writes one contig table per chain, for example
#   <prefix>_vdj_t_all_contig_annotations.csv
# contig_path is set in the config cell.
if (file.exists(contig_path)) {
  contigs <- read.csv(contig_path, stringsAsFactors = FALSE)
  # One row per cell: keep the productive chains and join their CDR3 calls.
  productive <- contigs[contigs$productive %in% c(TRUE, "True", "true"), ]
  clone_call <- tapply(
    productive$cdr3,
    productive$barcode,
    function(x) paste(sort(unique(x)), collapse = "_")
  )
  obj[["clone"]] <- unname(clone_call[colnames(obj)])

  # A clone that occupies many cells has expanded; sharing a clonotype across
  # clusters is direct evidence of a shared origin.
  clone_size <- table(obj[["clone", drop = TRUE]])
  print(sort(clone_size, decreasing = TRUE)[1:10])

  show_plot(
    DimPlot(obj, reduction = "umap", group.by = "clone") +
      ggplot2::theme(legend.position = "none"),
    "dimplot-clones", width = 7, height = 6
  )
} else {
  message("No contig file at '", contig_path, "'.")
  message("Set contig_path in the config cell to your Cell Ranger *_contig_annotations.csv.")
}