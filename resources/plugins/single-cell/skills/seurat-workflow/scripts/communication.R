#!/usr/bin/env Rscript
# Thin adapter for sc_communication: flags -> config globals -> fragments -> JSON report.
args <- commandArgs(trailingOnly = FALSE)
here <- sub("^--file=", "", args[startsWith(args, "--file=")])
here <- dirname(if (length(here)) here[[1]] else ".")
source(file.path(here, "_common.R"))
source(file.path(here, "_fragments.R"))

flags <- parse_flags()
input <- flag_req(flags, "input")
group_by <- flag_str(flags, "group_by", "cell_type")
interactions_path <- flag_str(flags, "interactions_path")
ccc_min_pct <- flag_num(flags, "min_pct", 0.10)
ccc_perms <- flag_int(flags, "n_perms", 200L)
ccc_fdr <- flag_num(flags, "fdr", 0.05)
ccc_output <- flag_str(flags, "interactions_output", "communication/interactions.tsv")
figures_dir <- flag_str(flags, "figures_dir")
assay <- flag_str(flags, "assay")
seed <- flag_int(flags, "seed", 42L)
min_genes <- NULL

data_path <- input
project_name <- flag_str(flags, "project", "phi")
output_object <- NULL

require_seurat()
tool_input <- input
tool_output <- NULL

sink(stderr())
run_fragments(c("setup", "load", "communication"))
sink()
source(file.path(here, "reports", "communication.R"))
