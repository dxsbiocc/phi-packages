# Group the embedding by the batch column. Cells of one donor sitting in one
# corner of the UMAP is the signature of an uncorrected batch effect.
batch_column <- sample_column

if (batch_column %in% colnames(obj[[]])) {
  DimPlot(obj, reduction = "umap", group.by = c("seurat_clusters", batch_column))
} else {
  message("No '", batch_column, "' column; add it to inspect the batch structure.")
}