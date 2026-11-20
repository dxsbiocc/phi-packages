#!/usr/bin/env Rscript
# Thin adapter for sc_annotate: flags -> config -> fragments -> JSON report.
args <- commandArgs(trailingOnly = FALSE)
here <- sub("^--file=", "", args[startsWith(args, "--file=")])
here <- dirname(if (length(here)) here[[1]] else ".")
source(file.path(here, "_common.R"))
source(file.path(here, "_fragments.R"))

flags <- parse_flags()
input <- flag_req(flags, "input")
output <- flag_str(flags, "output")
panel_path <- flag_str(flags, "panel")
mapping_path <- flag_str(flags, "mapping")

data_path <- input
project_name <- flag_str(flags, "project", "phi")
output_object <- output
figures_dir <- flag_str(flags, "figures_dir")
assay <- flag_str(flags, "assay")
seed <- flag_int(flags, "seed", 42L)
min_genes <- flag_num(flags, "feature_low")
group_by <- flag_str(flags, "group_by", "seurat_clusters")
label_column <- flag_str(flags, "label_column", "cell_type")
annotation_method <- flag_str(flags, "method", "panel")
min_score <- flag_num(flags, "min_score", 0)

require_seurat()
tool_input <- input
tool_output <- NULL

steps <- c("setup", "load")
if (!is.null(panel_path)) steps <- c(steps, "annotate")
if (!is.null(panel_path) || !is.null(mapping_path)) steps <- c(steps, "assign-labels")
if (!is.null(output)) steps <- c(steps, "save")

sink(stderr())
run_fragments(steps)
if (!is.null(output)) tool_output <- output
sink()
source(file.path(here, "reports", "annotate.R"))