#!/usr/bin/env python3
"""Thin adapter for sc_py_integrate: flags -> config globals -> fragments -> JSON report.

Guards that are notebook skips (missing batch column) are tool errors here: a
script run asked for integration explicitly, so silence would be a lie.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import emit, fail, flag_int, flag_nums, flag_req, flag_str, parse_flags
from _fragments import run_fragments, run_report

flags = parse_flags()
input_path = flag_req(flags, "input")
output = flag_str(flags, "output")
batch_key = flag_str(flags, "batch_key", "sample")
dims_use = flag_int(flags, "dims_use", 30)

namespace = {
    "data_path": input_path,
    "project_name": flag_str(flags, "project", "phi"),
    "output_object": output,
    "figures_dir": flag_str(flags, "figures_dir"),
    "min_genes": None,
    "batch_key": batch_key,
    "batch_method": flag_str(flags, "method", "harmony"),
    "n_pcs": dims_use,
    "dims": dims_use,
    "resolution": flag_nums(flags, "resolution", [0.8]),
    "run_umap": True,
    "run_tsne": False,
    "seed": flag_int(flags, "seed", 42),
    "plot_tag": "before",
    "emit": emit,
    "tool_input": input_path,
    "tool_output": None,
}

try:
    run_fragments(["setup", "load"], namespace)
    adata = namespace["adata"]
    if batch_key not in adata.obs.columns:
        fail(
            f"batch column '{batch_key}' not in obs; available: "
            + ", ".join(adata.obs.columns)
        )
    if "X_pca" not in adata.obsm:
        fail("no X_pca embedding; run sc_py_preprocess before sc_py_integrate")
    if adata.obs[batch_key].nunique() < 2:
        fail(f"batch column '{batch_key}' has a single level; nothing to integrate")
    run_fragments(["integrate", "batch-compare", "cluster"], namespace)
    if output:
        run_fragments(["save"], namespace)
        namespace["tool_output"] = output
except ValueError as error:
    fail(str(error))
run_report("integrate", namespace)