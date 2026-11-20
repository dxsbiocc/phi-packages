#!/usr/bin/env python3
"""Thin adapter for sc_py_preprocess: flags -> config globals -> fragments -> JSON report."""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import (
    emit,
    fail,
    flag_bool,
    flag_int,
    flag_num,
    flag_nums,
    flag_req,
    flag_str,
    flag_vec,
    parse_flags,
)
from _fragments import run_fragments, run_report

flags = parse_flags()
input_path = flag_req(flags, "input")
output = flag_str(flags, "output")
normalization = flag_str(flags, "normalization", "lognormalize")
if normalization != "lognormalize":
    fail(
        "--normalization " + normalization + " is not available on the scanpy route; "
        "SCTransform lives in the Seurat route (sc_preprocess)"
    )
n_pcs = flag_int(flags, "n_pcs", 30)

namespace = {
    "data_path": input_path,
    "project_name": flag_str(flags, "project", "phi"),
    "output_object": output,
    "figures_dir": flag_str(flags, "figures_dir"),
    "min_genes": None,
    "normalization": normalization,
    "scale_factor": flag_num(flags, "scale_factor", 10000),
    "n_hvg": flag_int(flags, "n_variable_features", 2000),
    "hvg_flavor": flag_str(flags, "hvg_flavor", "seurat_v3"),
    "regress_vars": flag_vec(flags, "regress_vars"),
    "n_pcs": n_pcs,
    "dims": flag_int(flags, "dims_use", n_pcs),
    "resolution": flag_nums(flags, "resolution", [0.8]),
    "run_umap": flag_bool(flags, "umap", True),
    "run_tsne": flag_bool(flags, "tsne", False),
    "batch_key": flag_str(flags, "batch_key"),
    "seed": flag_int(flags, "seed", 42),
    "plot_tag": "before",
    "emit": emit,
    "tool_input": input_path,
    "tool_output": None,
}

try:
    run_fragments(
        ["setup", "load", "normalize", "hvg", "scale", "pca", "cluster"], namespace
    )
    if output:
        run_fragments(["save"], namespace)
        namespace["tool_output"] = output
except ValueError as error:
    fail(str(error))
run_report("preprocess", namespace)