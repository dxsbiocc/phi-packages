#!/usr/bin/env Rscript
# Thin adapter for chip_annotate: annotate peaks against a user-supplied TxDb/GTF
# and draw genomic distribution + TSS coverage profiles. Accepts a peak file/dir,
# a sample sheet, or a saved .rds object (DBA or GRanges) from an earlier step.
args <- commandArgs(trailingOnly = FALSE)
here <- sub("^--file=", "", args[startsWith(args, "--file=")])
here <- dirname(if (length(here)) here[[1]] else ".")
source(file.path(here, "_common.R"))
source(file.path(here, "_fragments.R"))
source(file.path(here, "_config.R"))

require_chipseeker()
require_bioc("txdbmaker")
tool_input <- data_path
tool_output <- NULL

sink(stderr())
run_fragments("setup")
if (grepl("\\.rds$", data_path, ignore.case = TRUE)) {
  obj <- read_object(data_path)
  if (inherits(obj, "DBA")) {
    dba_object <- obj
    consensus_peaks <- DiffBind::dba.peakset(obj, bRetrieve = TRUE)
  } else if (inherits(obj, "GRanges")) {
    consensus_peaks <- obj
  } else if (inherits(obj, "GRangesList")) {
    peak_list <- obj
  } else {
    fail("saved object is not a DBA, GRanges or GRangesList; cannot annotate")
  }
} else {
  steps <- character(0)
  if (!is.null(sample_sheet_path) && nzchar(sample_sheet_path)) steps <- c(steps, "load-sample-sheet")
  steps <- c(steps, "load-peaks")
  run_fragments(steps)
}
steps <- c("annotate", "coverage-profile")
if (!is.null(output_object) && nzchar(output_object)) {
  steps <- c(steps, "save")
  tool_output <- output_object
}
run_fragments(steps)
sink()
source(file.path(here, "reports", "annotate.R"))