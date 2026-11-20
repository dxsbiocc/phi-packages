#!/usr/bin/env Rscript
# Report which single-cell R packages this environment provides, and which
# scNotebooks modules can therefore run here.
#
# scNotebooks installs its own packages inside Google Colab on every run
# (`install.packages`, `BiocManager::install`, `remotes::install_github`). Phi
# environments are read-only and prebuilt, so the honest first question is what
# the environment already has. This script answers it and maps the gaps to a
# concrete next step instead of failing mid-analysis.

args <- commandArgs(trailingOnly = FALSE)
here <- sub("^--file=", "", args[startsWith(args, "--file=")])
source(file.path(dirname(if (length(here)) here[[1]] else "."), "_common.R"))

flags <- parse_flags()
module <- flag_str(flags, "module")
write_panel <- flag_str(flags, "write_panel")

# package -> what it is needed for, and which scNotebooks modules rely on it.
catalogue <- list(
  list(package = "Seurat", role = "core single-cell object and workflow", modules = c(4, 5, 6, 7, 8, 9, 10, 11), core = TRUE),
  list(package = "SeuratObject", role = "Seurat class definitions", modules = c(4, 5, 6), core = TRUE),
  list(package = "SingleCellExperiment", role = "Bioconductor interchange container", modules = c(4, 7, 8), core = TRUE),
  list(package = "harmony", role = "batch correction (Module 05)", modules = c(4, 5, 6, 11), required = FALSE),
  list(package = "pheatmap", role = "marker and correlation heatmaps", modules = c(4, 10), required = FALSE),
  list(package = "clusterProfiler", role = "GO and pathway enrichment", modules = c(4, 10), required = FALSE),
  list(package = "org.Hs.eg.db", role = "human gene annotation for enrichment", modules = c(4, 10), required = FALSE),
  list(package = "ReactomePA", role = "Reactome pathway analysis (Module 10)", modules = c(10), required = FALSE),
  list(package = "enrichR", role = "Enrichr database queries (Module 10)", modules = c(10), required = FALSE),
  list(package = "SingleR", role = "reference-based cell type annotation", modules = c(4), required = FALSE),
  list(package = "celldex", role = "reference datasets for SingleR", modules = c(4), required = FALSE),
  list(package = "clustree", role = "cluster resolution tree", modules = c(10, 11), required = FALSE),
  list(package = "presto", role = "fast Wilcoxon marker test", modules = c(4, 10), required = FALSE),
  list(package = "miloR", role = "differential abundance (Milo)", modules = c(4), required = FALSE),
  list(package = "SoupX", role = "ambient RNA correction", modules = c(4), required = FALSE),
  list(package = "DropletUtils", role = "10x counts input and output for SoupX", modules = c(4), required = FALSE),
  list(package = "monocle3", role = "trajectory and pseudotime", modules = c(6), required = TRUE),
  list(package = "liana", role = "cell-cell communication", modules = c(7), required = TRUE),
  list(package = "scDblFinder", role = "doublet detection", modules = c(4), required = FALSE),
  list(package = "scater", role = "single-cell QC helpers", modules = c(4), required = FALSE),
  list(package = "hdf5r", role = "reading Cell Ranger .h5 matrices", modules = c(4, 8, 9, 10), required = FALSE),
  list(package = "semla", role = "spatial transcriptomics visualisation", modules = c(10), required = TRUE),
  list(package = "ArchR", role = "scATAC-seq analysis", modules = c(11), required = TRUE),
  list(package = "reticulate", role = "Python bridge used by ArchR and Module 12", modules = c(11, 12), required = FALSE),
  list(package = "uwot", role = "UMAP implementation used by Seurat (dependency)", modules = c(4, 5), core = TRUE),
  list(package = "Matrix", role = "sparse matrix support (Seurat dependency)", modules = c(4), core = TRUE)
)

probe <- function(entry) {
  installed <- requireNamespace(entry$package, quietly = TRUE)
  version <- if (installed) {
    tryCatch(as.character(utils::packageVersion(entry$package)), error = function(e) NA_character_)
  } else {
    NA_character_
  }
  list(
    package = entry$package,
    role = entry$role,
    installed = installed,
    version = version,
    core = isTRUE(entry$core),
    required = isTRUE(entry$required),
    modules = chr(as.integer(entry$modules))
  )
}

packages <- lapply(catalogue, probe)
names_tested <- vapply(packages, function(entry) entry$package, character(1))

# Module availability is stated as a fact about this environment, not a promise.
absent <- function(entries) {
  vapply(
    Filter(function(entry) !requireNamespace(entry$package, quietly = TRUE), entries),
    function(entry) entry$package,
    character(1)
  )
}

# A module's core path (the part these scripts implement) needs only its
# required packages; the rest are listed as optional sections that stay
# unavailable here. Saying "core_ready" avoids both overclaiming and refusing a
# module whose main workflow works.
modules <- lapply(1:13, function(number) {
  needed <- Filter(function(entry) number %in% entry$modules, catalogue)
  required <- Filter(function(entry) isTRUE(entry$required), needed)
  optional <- Filter(function(entry) !isTRUE(entry$required), needed)
  missing_required <- absent(required)
  list(
    module = number,
    packages_checked = length(needed),
    missing_required = chr(missing_required),
    missing_optional = chr(absent(optional)),
    core_ready = length(missing_required) == 0
  )
})

selected <- if (!is.null(module)) {
  Filter(function(entry) entry$module %in% as.integer(module), modules)
} else {
  modules
}

written <- NULL
if (!is.null(write_panel)) {
  # Emit a starter marker panel, so the annotate step can be tried immediately.
  # These are canonical symbols from the scNotebooks Module 04 manual panel.
  panel <- list(
    T_cell = c("CD3D", "CD3E", "TRAC", "TRBC1", "IL7R", "LTB"),
    B_cell = c("MS4A1", "CD79A", "CD79B", "CD74", "HLA-DRA", "CD19"),
    Plasma = c("MZB1", "XBP1", "JCHAIN", "SDC1", "TNFRSF17", "IGKC"),
    NK = c("NKG7", "KLRD1", "GNLY", "PRF1", "TRDC"),
    Myeloid = c("LYZ", "S100A8", "S100A9", "LGALS3", "CTSS", "MS4A7"),
    Dendritic = c("FCER1A", "CLEC10A", "CD1C", "ITGAX", "LILRA4", "IRF7"),
    Fibroblast = c("COL3A1", "PDGFRA", "LAMA2", "FBLN2", "DCN"),
    Epithelial = c("EPCAM", "KRT8", "KRT18", "KRT19", "CDH1"),
    Endothelial = c("PECAM1", "VWF", "CDH5", "CLDN5"),
    Erythrocyte = c("HBB", "HBA1", "HBA2", "ALAS2", "SLC4A1", "GYPA")
  )
  dir.create(dirname(write_panel), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(panel, write_panel, auto_unbox = FALSE, pretty = TRUE)
  info("panel ", write_panel)
  written <- write_panel
}

installed_count <- sum(vapply(packages, function(entry) entry$installed, logical(1)))

emit(list(
  environment = paste0("R ", getRversion()),
  packages_total = length(packages),
  packages_installed = installed_count,
  packages = packages,
  modules = selected,
  core_ready_modules = chr(vapply(
    Filter(function(entry) entry$core_ready, selected),
    function(entry) as.character(entry$module),
    character(1)
  )),
  blocked_modules = chr(vapply(
    Filter(function(entry) !entry$core_ready, selected),
    function(entry) as.character(entry$module),
    character(1)
  )),
  panel_written = written,
  notes = chr(c(
    if (installed_count < length(packages))
      "the plugin environment is read-only and prebuilt: a missing package needs a project environment (env_request), not an in-environment install",
    "scNotebooks installs packages inside Colab on every run; these scripts never do"
  ))
))
