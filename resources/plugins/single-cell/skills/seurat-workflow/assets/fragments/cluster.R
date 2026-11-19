# Graph-based clustering: link each cell to its nearest neighbours in PCA
# space, then find communities in that graph (Louvain by default; algorithm 4
# is Leiden, which guarantees well-connected communities).
obj <- FindNeighbors(obj, dims = seq_len(dims), verbose = FALSE)
analysis_steps <- c(analysis_steps, "neighbors")

# Resolution controls granularity: higher = more, smaller clusters. There is no
# single correct value; try several and judge against your biology. Each value
# gets its own metadata column so the choices stay comparable.
cluster_columns <- character(0)
for (value in resolution) {
  column <- if (length(resolution) == 1) "seurat_clusters" else paste0("seurat_clusters_", value)
  obj <- FindClusters(obj, resolution = value, algorithm = cluster_algorithm, verbose = FALSE)
  if (column != "seurat_clusters") {
    obj[[column]] <- obj[["seurat_clusters", drop = TRUE]]
  }
  cluster_columns <- c(cluster_columns, column)
}
analysis_steps <- c(analysis_steps, "clusters")

table(obj$seurat_clusters)

embeddings <- character(0)
if (run_umap) {
  obj <- RunUMAP(obj, dims = seq_len(dims), verbose = FALSE)
  embeddings <- c(embeddings, "umap")
}
if (run_tsne) {
  obj <- RunTSNE(obj, dims = seq_len(dims))
  embeddings <- c(embeddings, "tsne")
}
if (dims > n_pcs) {
  analysis_notes <- c(analysis_notes, sprintf("dims (%d) exceeds computed PCs (%d)", dims, n_pcs))
}

# A UMAP is a qualitative map: local neighbourhoods mean something, global
# distances and cluster sizes do not. Confirm any apparent structure with markers.
if (length(embeddings) > 0) {
  show_plot(DimPlot(obj, reduction = embeddings[[1]], label = TRUE) + NoLegend(),
            "dimplot-clusters", width = 7, height = 6)
}