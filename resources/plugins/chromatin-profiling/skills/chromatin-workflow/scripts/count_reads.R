#!/usr/bin/env Rscript
# Thin adapter for chip_count: count reads in a consensus peakset across samples
# and normalize for library size. --output saves the count matrix and consensus.
args <- commandArgs(trailingOnly = FALSE)
here <- sub("^--file=", "", args[startsWith(args, "--file=")])
here <- dirname(if (length(here)) here[[1]] else ".")
source(file.path(here, "_common.R"))
source(file.path(here, "_fragments.R"))
source(file.path(here, "_config.R"))

require_bioc("GenomicAlignments")
tool_input <- data_path
tool_output <- NULL

sink(stderr())
steps <- "setup"
if (!is.null(sample_sheet_path) && nzchar(sample_sheet_path)) steps <- c(steps, "load-sample-sheet")
steps <- c(steps, "load-peaks", "load-bam", "blacklist-filter", "consensus", "count", "normalize")
if (!is.null(output_object) && nzchar(output_object)) {
  steps <- c(steps, "save")
  tool_output <- output_object
}
run_fragments(steps)
sink()
source(file.path(here, "reports", "count.R"))