#!/usr/bin/env Rscript
# Thin adapter for sc_preprocess: flags -> config -> fragments -> JSON report.
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
assay <- flag_str(flags, "assay")
seed <- flag_int(flags, "seed", 42L)
min_genes <- flag_num(flags, "feature_low")
normalization <- flag_str(flags, "normalization", "lognormalize")
scale_factor <- flag_num(flags, "scale_factor", 10000)
selection_method <- flag_str(flags, "selection_method", "vst")
n_hvg <- flag_int(flags, "n_variable_features", 2000L)
n_pcs <- flag_int(flags, "n_pcs", 30L)
dims <- flag_int(flags, "dims_use", n_pcs)
resolution <- flag_nums(flags, "resolution", 0.8)
cluster_algorithm <- flag_int(flags, "cluster_algorithm", 1L)
run_umap <- flag_bool(flags, "umap", TRUE)
run_tsne <- flag_bool(flags, "tsne", FALSE)
regress_vars <- flag_vec(flags, "regress_vars")

require_seurat()
tool_input <- input
tool_output <- NULL

sink(stderr())
run_fragments(c("setup", "load", "normalize", "hvg", "scale", "pca", "cluster"))
if (!is.null(output)) {
  run_fragments("save")
  tool_output <- output
}
sink()
source(file.path(here, "reports", "preprocess.R"))