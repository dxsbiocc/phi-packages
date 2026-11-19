# Seurat's VlnPlot(ncol = 3) combines panels with patchwork's `&` operator,
# which raises "Can't find method for generic `&`" on the ggplot2 4.x that
# this environment ships. Building one plot per metric and assembling them with wrap_plots
# avoids that operator entirely.
qc_panels <- function(object, features, ncol = 3) {
  features <- features[features %in% colnames(object[[]])]
  panels <- lapply(features, function(feature) {
    VlnPlot(object, features = feature, pt.size = 0) +
      ggplot2::theme(legend.position = "none") +
      ggplot2::labs(x = NULL, title = feature)
  })
  patchwork::wrap_plots(panels, ncol = ncol)
}

# plot_tag names the pass ("before"/"after"): the script tools run this fragment
# twice, once on each side of the filter, so the two violins are comparable.
show_plot(qc_panels(obj, c("nFeature_RNA", "nCount_RNA", "percent.mt")),
          paste0("qc-violin-", plot_tag), width = 11, height = 4)

# Joint views catch what one-dimensional plots hide: a cell with many genes and
# a high mitochondrial fraction is stressed, not healthy-and-busy.
show_plot(FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mt"),
          paste0("qc-scatter-mt-", plot_tag))

show_plot(FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA"),
          paste0("qc-scatter-genes-", plot_tag))