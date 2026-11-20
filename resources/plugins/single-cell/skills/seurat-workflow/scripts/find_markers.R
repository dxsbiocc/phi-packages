#!/usr/bin/env Rscript
# Thin adapter for sc_markers: flags -> config -> fragments -> JSON report.
args <- commandArgs(trailingOnly = FALSE)
here <- sub("^--file=", "", args[startsWith(args, "--file=")])
here <- dirname(if (length(here)) here[[1]] else ".")
source(file.path(here, "_common.R"))
source(file.path(here, "_fragments.R"))

flags <- parse_flags()
input <- flag_req(flags, "input")
output <- flag_str(flags, "output")

data_path <- input
project_name <- flag_str(flags, "project", "phi")
output_object <- output
figures_dir <- flag_str(flags, "figures_dir")
tables_dir <- flag_str(flags, "tables_dir")
assay <- flag_str(flags, "assay")
seed <- flag_int(flags, "seed", 42L)
min_genes <- flag_num(flags, "feature_low")
group_by <- flag_str(flags, "group_by", "seurat_clusters")
test_use <- flag_str(flags, "test", "wilcox")
logfc_threshold <- flag_num(flags, "logfc_threshold", 0.25)
min_pct <- flag_num(flags, "min_pct", 0.1)
only_positive <- flag_bool(flags, "only_positive", TRUE)
top_n <- flag_int(flags, "top_n", 10L)
max_cells_per_ident <- flag_int(flags, "max_cells_per_ident")
plot_heatmap <- flag_bool(flags, "heatmap", TRUE)

require_seurat()
tool_input <- input
tool_output <- NULL

sink(stderr())
run_fragments(c("setup", "load", "markers"))
if (!is.null(output)) {
  run_fragments("save")
  tool_output <- output
}
sink()
source(file.path(here, "reports", "markers.R"))