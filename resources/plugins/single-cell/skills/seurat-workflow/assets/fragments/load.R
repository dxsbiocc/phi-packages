# data_path and project_name are set in the config cell above.
# A NULL min_genes means "no creation gate" (the script tools pass it that way).
gate <- if (is.null(min_genes)) 0 else min_genes
if (dir.exists(data_path)) {
  # Cell Ranger writes features.tsv.gz + barcodes.tsv.gz + matrix.mtx.gz here.
  counts <- Read10X(data.dir = data_path)

  # HARD GATE: min.features drops cells, and if it drops every cell Seurat fails
  # later with the opaque "LogMap object: Rownames must be supplied". Check the
  # distribution first, then filter at object creation.
  detected <- Matrix::colSums(counts > 0)
  cat(sprintf("genes detected per cell: min %d, median %d, max %d\n",
              min(detected), as.integer(median(detected)), max(detected)))

  if (sum(detected > gate) == 0) {
    stop(sprintf(paste(
      "No cell has more than %d detected genes, so min.features would remove",
      "every cell. This usually means a filtered/subset matrix, a different",
      "feature ID type, or an unusually shallow library. Lower min_genes in the",
      "config cell or check that the matrix rows are genes."
    ), gate), call. = FALSE)
  }

  obj <- CreateSeuratObject(
    counts = counts,
    project = project_name,
    min.cells = 3,          # drop genes seen in fewer than 3 cells
    min.features = gate  # drop cells with fewer genes than the config value
  )
  rm(counts)
} else if (grepl("[.]rds$", data_path)) {
  obj <- readRDS(data_path)
  if (!inherits(obj, "Seurat")) {
    stop("The .rds does not contain a Seurat object: ", paste(class(obj), collapse = "/"))
  }
} else if (grepl("[.]csv$|[.]tsv$", data_path)) {
  # A plain counts table: genes as rows, cells as columns, first column gene id.
  frame <- read.delim(data_path, sep = if (grepl("[.]csv$", data_path)) "," else "\t",
                      row.names = 1, check.names = FALSE)
  obj <- CreateSeuratObject(
    counts = as.matrix(frame),
    project = project_name,
    min.cells = 3,
    min.features = gate
  )
  rm(frame)
} else if (grepl("[.]h5$", data_path)) {
  stop(paste(
    "A Cell Ranger .h5 needs the hdf5r package, which the plugin environment does not ship.",
    "Pass the matching filtered_feature_bc_matrix directory instead."
  ), call. = FALSE)
} else {
  stop("Unsupported input: ", data_path, call. = FALSE)
}

# An explicit assay choice applies to every step below (markers, pseudobulk...).
if (!is.null(assay) && nzchar(assay)) {
  if (!assay %in% names(obj@assays)) {
    stop("assay not found: ", assay, "; available: ", paste(names(obj@assays), collapse = ", "))
  }
  DefaultAssay(obj) <- assay
}

# Genes (features) x cells is the Seurat convention.
loaded_assay <- DefaultAssay(obj)
cat(sprintf("%d genes x %d cells\n", nrow(obj), ncol(obj)))
cat("Assays:", paste(names(obj@assays), collapse = ", "), "\n")
cat("Layers:", paste(Layers(obj[[DefaultAssay(obj)]]), collapse = ", "), "\n")