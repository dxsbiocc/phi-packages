# Raw UMI counts are not comparable across cells: a cell with twice the reads
# shows twice the counts for every gene, which is a technical artifact, not
# biology. LogNormalize divides each cell by its total, rescales to scale_factor,
# then log-transforms to stabilise variance.
if (identical(normalization, "sct")) {
  # SCTransform replaces NormalizeData + FindVariableFeatures + ScaleData with
  # one regularised negative-binomial regression per gene; it models sequencing
  # depth instead of dividing it out, and adds an "SCT" assay whose scale.data
  # holds Pearson residuals. The hvg and scale cells below then skip themselves.
  layers <- tryCatch(Layers(obj[[DefaultAssay(obj)]]), error = function(e) character(0))
  if (!"counts" %in% layers) {
    analysis_notes <- c(analysis_notes, paste0(
      "assay '", DefaultAssay(obj), "' has no counts layer (found: ",
      paste(layers, collapse = ", "), "); normalization assumes raw counts"
    ))
  }
  obj <- SCTransform(obj, verbose = FALSE, assay = DefaultAssay(obj))
  analysis_steps <- c(analysis_steps, "sctransform")
} else {
  layers <- tryCatch(Layers(obj[[DefaultAssay(obj)]]), error = function(e) character(0))
  if (!"counts" %in% layers) {
    analysis_notes <- c(analysis_notes, paste0(
      "assay '", DefaultAssay(obj), "' has no counts layer (found: ",
      paste(layers, collapse = ", "), "); normalization assumes raw counts"
    ))
  }
  obj <- NormalizeData(
    obj,
    normalization.method = "LogNormalize",
    scale.factor = scale_factor,
    verbose = FALSE
  )
  analysis_steps <- c(analysis_steps, "normalize")
}

obj[[DefaultAssay(obj)]]$data[1:5, 1:3]