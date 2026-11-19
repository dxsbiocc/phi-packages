#!/usr/bin/env Rscript
# Inspect a single-cell input and report what it contains.
#
# Accepts a Seurat .rds, a 10x Cell Ranger matrix directory
# (filtered_feature_bc_matrix or raw_feature_bc_matrix), a .mtx matrix, or a
# dense .csv/.tsv count table. Writes a Seurat object when --output is given.
#
# Ported from scNotebooks Module 04 "Data input, sparse UMI-count matrix"
# and Module 09 "Read in raw gene counts and metadata".

args <- commandArgs(trailingOnly = FALSE)
here <- sub("^--file=", "", args[startsWith(args, "--file=")])
source(file.path(dirname(if (length(here)) here[[1]] else "."), "_common.R"))

flags <- parse_flags()
input <- flag_req(flags, "input")
output <- flag_str(flags, "output")
project <- flag_str(flags, "project", basename(normalizePath(input, mustWork = FALSE)))
min_cells <- flag_int(flags, "min_cells", 0L)
min_features <- flag_int(flags, "min_features", 0L)
var_names <- flag_str(flags, "var_names", "gene_symbols")

require_seurat()

detect_kind <- function(path) {
  if (dir.exists(path)) {
    markers <- c("features.tsv.gz", "genes.tsv.gz", "features.tsv")
    has_markers <- any(file.exists(file.path(path, markers)))
    has_matrix <- file.exists(file.path(path, "matrix.mtx.gz")) ||
      file.exists(file.path(path, "matrix.mtx"))
    if (has_markers && has_matrix) return("10x-matrix-dir")
    # A Cell Ranger `outs/` directory: look one level down.
    nested <- file.path(path, "filtered_feature_bc_matrix")
    if (dir.exists(nested)) return("10x-matrix-dir")
    return("unknown-directory")
  }
  lower <- tolower(path)
  if (endsWith(lower, ".rds")) return("seurat-rds")
  if (endsWith(lower, ".h5")) return("10x-h5")
  if (endsWith(lower, ".h5ad")) return("anndata")
  if (endsWith(lower, ".mtx") || endsWith(lower, ".mtx.gz")) return("matrix-market")
  if (endsWith(lower, ".csv") || endsWith(lower, ".tsv") || endsWith(lower, ".txt")) return("dense-table")
  "unknown-file"
}

load_input <- function(path, kind) {
  switch(
    kind,
    "seurat-rds" = read_object(path),
    "10x-matrix-dir" = {
      directory <- path
      nested <- file.path(path, "filtered_feature_bc_matrix")
      if (dir.exists(nested)) directory <- nested
      info("reading 10x matrix directory ", directory)
      Seurat::CreateSeuratObject(
        counts = Seurat::Read10X(data.dir = directory, gene.column = 1),
        project = project,
        min.cells = min_cells,
        min.features = min_features
      )
    },
    "matrix-market" = {
      info("reading matrix market ", path)
      counts <- utils::read.delim(path, header = FALSE, comment.char = "%")
      fail(
        "--input ", path, " is a bare Matrix Market file without feature/barcode ",
        "names; pass its 10x directory instead so cell and gene identifiers survive"
      )
    },
    "dense-table" = {
      info("reading dense table ", path)
      separator <- if (endsWith(tolower(path), ".csv")) "," else "\t"
      frame <- utils::read.csv(path, sep = separator, row.names = 1, check.names = FALSE)
      counts <- as.matrix(frame)
      storage.mode(counts) <- "numeric"
      Seurat::CreateSeuratObject(
        counts = Matrix::Matrix(counts, sparse = TRUE),
        project = project,
        min.cells = min_cells,
        min.features = min_features
      )
    },
    "10x-h5" = fail(
      "--input is a 10x .h5 file; reading it needs the hdf5r package, which the ",
      "plugin environment does not ship. Convert it first (Cell Ranger writes the ",
      "same matrix as an mtx directory) or pass that directory with --input"
    ),
    "anndata" = fail(
      "--input is .h5ad; Seurat cannot read AnnData. Use the scanpy skill, or ",
      "convert to .rds with the cross-route notes in scanpy-workflow/references/cross-route-handoff.md"
    ),
    "unknown-directory" = fail(
      "--input is a directory with no 10x matrix files (matrix.mtx + features.tsv); ",
      "expected a Cell Ranger filtered_feature_bc_matrix or raw_feature_bc_matrix folder"
    ),
    "unknown-file" = fail(
      "unsupported --input; expected a Seurat .rds, a 10x matrix directory, ",
      "or a .csv/.tsv count table"
    )
  )
}

kind <- detect_kind(input)
info("input kind: ", kind)
object <- load_input(input, kind)

assays <- names(object@assays)
default_assay <- Seurat::DefaultAssay(object)
layers <- tryCatch(SeuratObject::Layers(object[[default_assay]]), error = function(e) character(0))
metadata_columns <- colnames(object[[]])
reductions <- SeuratObject::Reductions(object)
graphs <- tryCatch(SeuratObject::Graphs(object), error = function(e) character(0))

# Counts kept alongside the normalized data are what makes a downstream
# differential-expression handoff possible, so report them explicitly.
matrix_shape <- list(
  genes = nrow(object),
  cells = ncol(object),
  assays = chr(assays),
  default_assay = default_assay,
  layers = chr(layers),
  reductions = chr(reductions),
  graphs = chr(graphs)
)

has_counts <- "counts" %in% layers
value_class <- tryCatch(
  class(SeuratObject::LayerData(object[[default_assay]], layer = layers[[1]]))[[1]],
  error = function(e) NA_character_
)

# Metadata columns worth handing to a downstream tool, with their spread.
interesting <- intersect(
  metadata_columns,
  c(
    "orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", "percent.hb",
    "seurat_clusters", "seurat_annotations", "cell_type", "celltype",
    "sample", "condition", "group", "stim", "batch"
  )
)
metadata <- lapply(interesting, function(column) {
  values <- object[[column, drop = TRUE]]
  if (is.numeric(values)) {
    list(
      column = column,
      kind = "number",
      min = round(min(values, na.rm = TRUE), 4),
      max = round(max(values, na.rm = TRUE), 4),
      median = round(stats::median(values, na.rm = TRUE), 4),
      missing = sum(is.na(values))
    )
  } else {
    categories <- sorted_levels(values)
    list(
      column = column,
      kind = "category",
      distinct = length(categories),
      levels = chr(utils::head(categories, 40)),
      missing = sum(is.na(values))
    )
  }
})

ready_for <- character(0)
if (has_counts) ready_for <- c(ready_for, "qc")
if (length(reductions) > 0) ready_for <- c(ready_for, "markers")
if ("seurat_clusters" %in% metadata_columns) ready_for <- c(ready_for, "annotate")
if (has_counts && length(metadata_columns) > 0) ready_for <- c(ready_for, "pseudobulk")

written <- NULL
if (!is.null(output)) {
  write_object(object, output)
  written <- output
}

emit(list(
  input = input,
  kind = kind,
  output = written,
  shape = matrix_shape,
  counts_present = has_counts,
  value_class = value_class,
  metadata = metadata,
  metadata_columns = chr(metadata_columns),
  candidate_steps = chr(ready_for),
  notes = chr(c(
    if (!has_counts) "no counts layer: re-import the raw matrix before QC or pseudobulk",
    if (!is.na(value_class) && value_class == "dgCMatrix") "sparse matrix retained",
    if (kind == "10x-matrix-dir") "Cell Ranger matrix read with gene symbols from features.tsv"
  ))
))
