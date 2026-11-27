# Packages this notebook uses. The plugin environment already provides Seurat 6
# and SingleCellExperiment; the optional ones are checked up front so a gap shows
# here rather than halfway through the analysis.
have <- function(pkg) requireNamespace(pkg, quietly = TRUE)

need <- function(...) {
  for (pkg in c(...)) {
    if (!have(pkg)) {
      stop(sprintf("Package '%s' is not installed in this environment.", pkg),
           call. = FALSE)
    }
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
}

need("Seurat")

# Large objects exceed Seurat's default 500 MB future limit.
options(future.globals.maxSize = 2 * 1024^3)
set.seed(seed)

# Figures render inline; the size keeps verbose legends readable.
options(repr.plot.width = 8, repr.plot.height = 5)

# Bookkeeping shared with the script tools: the same fragments run there, and
# the report layer reads these variables to build the JSON result.
figure_paths <- list()
analysis_steps <- character(0)
analysis_notes <- character(0)
qc_scored <- NULL
qc_filters <- list()
obj_before <- NULL
doublet_result <- NULL
pca_variance <- NULL
cluster_columns <- character(0)
embeddings <- list()
markers <- NULL
marker_tables <- list()
annotation_panel <- NULL
annotation_scores <- NULL
annotation_assignments <- NULL
annotation_automated <- NULL
pb_info <- NULL
ccc_info <- NULL

# Show a plot inline and, when figures_dir is set, also write a PNG so a run
# leaves reviewable figures behind. Returns the path (or NULL) invisibly.
show_plot <- function(plot, name, width = 8, height = 5) {
  # Inline display only inside a kernel: a bare print() under Rscript would
  # open the default PDF device and leave Rplots.pdf in the working directory.
  if ("IRkernel" %in% loadedNamespaces()) print(plot)
  path <- NULL
  if (!is.null(figures_dir) && nzchar(figures_dir)) {
    dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
    path <- file.path(figures_dir, paste0(name, ".png"))
    device <- if (have("ragg")) ragg::agg_png else grDevices::png
    device(path, width = width, height = height, units = "in", res = 150)
    print(plot)
    grDevices::dev.off()
    figure_paths[[name]] <<- path
  }
  invisible(path)
}

cat("Seurat", as.character(packageVersion("Seurat")), "on R", as.character(getRversion()), "\n")