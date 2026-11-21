#!/usr/bin/env Rscript
# Thin adapter for chip_qc: FRiP, replicate reproducibility and blacklist filtering.
# CLI flags become the notebook config variables, the fragment library runs the
# analysis, reports/qc.R emits the JSON contract.
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
steps <- c(steps, "load-peaks", "load-bam", "qc-frip", "qc-reproducibility", "blacklist-filter")
run_fragments(steps)
sink()
source(file.path(here, "reports", "qc.R"))