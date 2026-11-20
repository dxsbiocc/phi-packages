# Only genes that vary across cells carry information about cell identity.
# Housekeeping genes vary little and would add noise to the distance metric.
# vst fits a mean-variance trend on the counts and ranks genes by standardized
# variance, so highly expressed genes do not win just for being loud.
if (DefaultAssay(obj) == "SCT") {
  message("SCT assay already carries variable features; skipping FindVariableFeatures.")
} else {
  obj <- FindVariableFeatures(
    obj,
    selection.method = selection_method,
    nfeatures = n_hvg,
    verbose = FALSE
  )
  analysis_steps <- c(analysis_steps, "variable-features")
}

# The 10 most variable genes are usually ribosomal/mitochondrial contamination
# or a single dominant cell type. Look before trusting them.
head(VariableFeatures(obj), 10)

show_plot(VariableFeaturePlot(obj), "hvg-mean-variance")