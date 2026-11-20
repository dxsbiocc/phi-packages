# PCA finds the directions of greatest variance. The first components usually
# separate major cell types; later ones carry noise — but rare populations can
# hide in later PCs, so keep a few more than the elbow suggests.
obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = n_pcs, verbose = FALSE)
analysis_steps <- c(analysis_steps, "pca")

show_plot(DimPlot(obj, reduction = "pca"), "pca-dimplot")

# Variance explained: the elbow tells you how many components are worth keeping.
# Seurat keeps only the computed PCs, so shares are measured against those.
pca_stdev <- Stdev(obj, reduction = "pca")
pca_share <- 100 * (pca_stdev^2) / sum(pca_stdev^2)
pca_variance <- list(
  pcs = length(pca_stdev),
  variance_explained = lapply(seq_along(pca_stdev), function(i) round(pca_share[[i]], 3)),
  first_ten_percent = round(sum(pca_share[seq_len(min(10, length(pca_share)))]), 3)
)
round(pca_share, 2)

show_plot(ElbowPlot(obj, ndims = length(pca_stdev)), "pca-elbow", width = 7, height = 5)