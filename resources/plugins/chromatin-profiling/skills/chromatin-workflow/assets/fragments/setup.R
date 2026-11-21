# Packages this notebook uses. The plugin environment already provides the
# ChIP/CUT&Tag downstream stack (DiffBind, ChIPseeker, csaw, DESeq2/edgeR and the
# GenomicRanges IO core); they are attached up front so a gap shows here rather
# than halfway through the analysis. Nothing is installed: the environment is
# read-only and prebuilt.
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

need("GenomicRanges", "rtracklayer", "DiffBind", "ChIPseeker",
     "GenomicAlignments", "Rsamtools", "GenomicFeatures",
     "DESeq2", "edgeR", "ComplexHeatmap", "ggplot2")

# csaw (window-based differential binding for broad marks) has no osx-arm64
# bioconda build, so it is optional: present on linux-x64 / darwin-x64, a reported
# gap on Apple Silicon. The csaw fragment skips gracefully when it is absent.
if (!have("csaw")) {
  cat("note: csaw is not available on this platform (no osx-arm64 build); ",
      "broad-mark window analysis is a gap, DiffBind still works\n", sep = "")
}

# DiffBind's counting and csaw's windowCounts parallelise over BAM files; cap the
# workers so a notebook does not fork the whole machine.
options(mc.cores = 2)
set.seed(seed)

# Figures render inline; the size keeps verbose legends readable.
options(repr.plot.width = 8, repr.plot.height = 5)

# Under Rscript (script tools and --self-check) some package calls (e.g.
# DiffBind::dba.plotMA) open the default device as a side effect and litter the
# working directory with Rplots.pdf. Route such incidental draws to a null device;
# show_plot still opens its own PNG (or the kernel display) for the real figures.
if (!"IRkernel" %in% loadedNamespaces()) grDevices::pdf(NULL)

# Bookkeeping shared with the script tools: the same fragments run there, and
# the report layer reads these variables to build the JSON result.
figure_paths <- list()
analysis_steps <- character(0)
analysis_notes <- character(0)
peak_list <- NULL
consensus_peaks <- NULL
occupancy_matrix <- NULL
count_matrix <- NULL
normalized_counts <- NULL
dba_object <- NULL
csaw_result <- NULL
diff_results <- NULL
qc_metrics <- NULL
blacklist_removed <- NULL
annotation_result <- NULL
anno_df <- NULL
txdb <- NULL
sample_sheet <- NULL
pca_result <- NULL

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

cat("DiffBind", as.character(packageVersion("DiffBind")),
    "| ChIPseeker", as.character(packageVersion("ChIPseeker")),
    if (have("csaw")) paste("| csaw", as.character(packageVersion("csaw"))) else "| csaw absent",
    "on R", as.character(getRversion()), "\n")