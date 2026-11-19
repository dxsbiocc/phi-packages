#!/usr/bin/env python3
"""Thin adapter for sc_py_pseudobulk: flags -> config globals -> fragments -> JSON report."""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import emit, fail, flag_int, flag_num, flag_req, flag_str, parse_flags
from _fragments import run_fragments, run_report

flags = parse_flags()
input_path = flag_req(flags, "input")
sample_column = flag_req(flags, "sample_column")

namespace = {
    "data_path": input_path,
    "project_name": flag_str(flags, "project", "phi"),
    "min_genes": None,
    "sample_column": sample_column,
    "group_column": flag_str(flags, "group_column"),
    "min_cells": flag_int(flags, "min_cells", 10),
    "min_samples_per_group": flag_int(flags, "min_samples_per_group", 2),
    "min_count": flag_num(flags, "min_count", 1),
    "min_cells_expressing": flag_int(flags, "min_cells_expressing", 3),
    "counts_output": flag_str(flags, "counts_output"),
    "samples_output": flag_str(flags, "samples_output"),
    "seed": flag_int(flags, "seed", 42),
    "emit": emit,
    "tool_input": input_path,
    "tool_output": None,
}

try:
    run_fragments(["setup", "load"], namespace)
    if sample_column not in namespace["adata"].obs.columns:
        fail(
            f"sample column '{sample_column}' not in obs; available: "
            + ", ".join(namespace["adata"].obs.columns)
        )
    run_fragments(["pseudobulk"], namespace)
except ValueError as error:
    fail(str(error))
run_report("pseudobulk", namespace)