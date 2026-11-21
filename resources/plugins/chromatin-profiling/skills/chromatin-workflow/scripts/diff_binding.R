#!/usr/bin/env Rscript
# Thin adapter for chip_diff: differential binding across conditions.
# --method deseq2|edger uses DiffBind on a consensus peakset; --method csaw uses
# csaw window-based analysis for broad marks (csaw is a gap on darwin-arm64).
args <- commandArgs(trailingOnly = FALSE)
here <- sub("^--file=", "", args[startsWith(args, "--file=")])
here <- dirname(if (length(here)) here[[1]] else ".")
source(file.path(here, "_common.R"))
source(file.path(here, "_fragments.R"))
source(file.path(here, "_config.R"))

method <- flag_str(flags, "method", "deseq2")
if (!method %in% c("deseq2", "edger", "csaw")) {
  fail(paste0("unknown --method '", method, "'; expected deseq2, edger or csaw"))
}
if (method == "csaw") require_csaw() else require_diffbind()
tool_input <- data_path
tool_output <- NULL

sink(stderr())
if (method == "csaw") {
  steps <- c("setup", "load-sample-sheet", "load-bam", "csaw")
} else {
  steps <- c(
    "setup", "load-sample-sheet", "load-peaks", "load-bam",
    "blacklist-filter", "consensus", "count",
    if (method == "edger") "diff-edger" else "diff-deseq2"
  )
}
if (!is.null(output_object) && nzchar(output_object)) {
  steps <- c(steps, "save")
  tool_output <- output_object
}
run_fragments(steps)
sink()
source(file.path(here, "reports", "diff.R"))