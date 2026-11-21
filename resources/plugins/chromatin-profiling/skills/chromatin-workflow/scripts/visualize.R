#!/usr/bin/env Rscript
# Thin adapter for chip_visualize: exploratory views of a peak x sample count
# matrix — sample PCA and a peak-overlap UpSet plot. Entry point is either a count
# table (--counts or a .csv/.tsv input) or peaks + BAMs that get counted here.
args <- commandArgs(trailingOnly = FALSE)
here <- sub("^--file=", "", args[startsWith(args, "--file=")])
here <- dirname(if (length(here)) here[[1]] else ".")
source(file.path(here, "_common.R"))
source(file.path(here, "_fragments.R"))
source(file.path(here, "_config.R"))

require_bioc("ggplot2")
tool_input <- data_path
tool_output <- NULL

sink(stderr())
steps <- "setup"
count_input <- (!is.null(count_path) && nzchar(count_path)) ||
  (is.null(sample_sheet_path) &&
    grepl("\\.(csv|tsv|txt)$", data_path, ignore.case = TRUE) &&
    file.exists(data_path))
if (count_input) {
  steps <- c(steps, "load-counts")
} else {
  if (!is.null(sample_sheet_path) && nzchar(sample_sheet_path)) steps <- c(steps, "load-sample-sheet")
  steps <- c(steps, "load-peaks", "load-bam", "blacklist-filter", "consensus", "count")
}
steps <- c(steps, "normalize", "pca", "upset")
if (!is.null(output_object) && nzchar(output_object)) {
  steps <- c(steps, "save")
  tool_output <- output_object
}
run_fragments(steps)
sink()
source(file.path(here, "reports", "visualize.R"))