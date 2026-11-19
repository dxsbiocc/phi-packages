#!/usr/bin/env Rscript
# Thin adapter for sc_qc: CLI flags become the notebook config variables, the
# fragment library runs the analysis, reports/qc.R emits the JSON contract.
args <- commandArgs(trailingOnly = FALSE)
here <- sub("^--file=", "", args[startsWith(args, "--file=")])
here <- dirname(if (length(here)) here[[1]] else ".")
source(file.path(here, "_common.R"))
source(file.path(here, "_fragments.R"))

flags <- parse_flags()
input <- flag_req(flags, "input")
output <- flag_str(flags, "output")
save_object <- flag_bool(flags, "save_object", TRUE)

data_path <- input
project_name <- flag_str(flags, "project", "phi")
output_object <- output
figures_dir <- flag_str(flags, "figures_dir")
plot_tag <- "before"
assay <- flag_str(flags, "assay")
seed <- flag_int(flags, "seed", 42L)
min_genes <- flag_num(flags, "feature_low")
max_genes <- flag_num(flags, "feature_high")
min_counts <- flag_num(flags, "count_low")
max_counts <- flag_num(flags, "count_high")
mt_max <- flag_num(flags, "mt_max", 10)
mt_pattern <- flag_str(flags, "mt_pattern", "^MT-")
ribo_pattern <- flag_str(flags, "ribo_pattern", "^RP[SL]")
hb_pattern <- flag_str(flags, "hb_pattern", "^HB[ABDEGMQZ][0-9]*$")
quantiles <- flag_nums(flags, "quantiles", c(0.01, 0.99))
sample_column <- flag_str(flags, "sample_column", "orig.ident")
detect_doublets <- flag_bool(flags, "detect_doublets", FALSE)
doublet_expected_rate <- flag_num(flags, "doublet_expected_rate")

require_seurat()
tool_input <- input
tool_output <- NULL

sink(stderr())
run_fragments(c("setup", "load", "qc-metrics", "doublets", "qc-plots"))
plot_tag <- "after"
run_fragments(c("filter", "qc-plots"))
if (!is.null(output) && save_object) {
  run_fragments("save")
  tool_output <- output
}
sink()
source(file.path(here, "reports", "qc.R"))