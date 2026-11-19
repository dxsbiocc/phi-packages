#!/usr/bin/env Rscript
# Thin adapter for sc_pseudobulk: flags -> config -> fragments -> JSON report.
args <- commandArgs(trailingOnly = FALSE)
here <- sub("^--file=", "", args[startsWith(args, "--file=")])
here <- dirname(if (length(here)) here[[1]] else ".")
source(file.path(here, "_common.R"))
source(file.path(here, "_fragments.R"))

flags <- parse_flags()
input <- flag_req(flags, "input")
sample_column <- flag_req(flags, "sample_column")
group_column <- flag_str(flags, "group_column")
counts_output <- flag_str(flags, "counts_output")
samples_output <- flag_str(flags, "samples_output")
min_cells <- flag_int(flags, "min_cells", 10L)
min_samples_per_group <- flag_int(flags, "min_samples_per_group", 2L)
min_count <- flag_num(flags, "min_count", 1)
min_cells_expressing <- flag_int(flags, "min_cells_expressing", 3L)

data_path <- input
project_name <- flag_str(flags, "project", "phi")
output_object <- NULL
figures_dir <- NULL
assay <- flag_str(flags, "assay")
seed <- flag_int(flags, "seed", 42L)
min_genes <- flag_num(flags, "feature_low")

require_seurat()
tool_input <- input
tool_output <- NULL

sink(stderr())
run_fragments(c("setup", "load"))
resolve_group(obj, sample_column)
run_fragments("pseudobulk")
sink()
source(file.path(here, "reports", "pseudobulk.R"))