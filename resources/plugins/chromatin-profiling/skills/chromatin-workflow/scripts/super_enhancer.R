#!/usr/bin/env Rscript
# Thin adapter for chip_super_enhancer: flags -> config -> fragments -> JSON report.
args <- commandArgs(trailingOnly = FALSE)
here <- sub("^--file=", "", args[startsWith(args, "--file=")])
here <- dirname(if (length(here)) here[[1]] else ".")
source(file.path(here, "_common.R"))
source(file.path(here, "_fragments.R"))
source(file.path(here, "_config.R"))

require_bioc("GenomicRanges")
tool_input <- data_path
tool_output <- NULL

sink(stderr())
run_fragments(c("setup", "load-peaks", "super-enhancer"))
sink()
source(file.path(here, "reports", "super_enhancer.R"))
