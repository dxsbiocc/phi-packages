#!/usr/bin/env python3
"""Thin adapter for sc_py_communication: flags -> config globals -> fragments -> JSON report."""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import emit, fail, flag_int, flag_num, flag_req, flag_str, parse_flags
from _fragments import run_fragments, run_report

flags = parse_flags()
input_path = flag_req(flags, "input")
output = flag_str(flags, "output")

namespace = {
    "data_path": input_path,
    "project_name": flag_str(flags, "project", "phi"),
    "output_object": output,
    "figures_dir": flag_str(flags, "figures_dir"),
    "min_genes": None,
    "group_by": flag_str(flags, "group_by", "cell_type"),
    "interactions_path": flag_str(flags, "interactions_path"),
    "ccc_min_pct": flag_num(flags, "min_pct", 0.10),
    "ccc_perms": flag_int(flags, "n_perms", 500),
    "ccc_fdr": flag_num(flags, "fdr", 0.05),
    "ccc_output": flag_str(flags, "interactions_output"),
    "seed": flag_int(flags, "seed", 42),
    "plot_tag": "before",
    "emit": emit,
    "tool_input": input_path,
    "tool_output": None,
}

try:
    run_fragments(["setup", "load", "communication"], namespace)
    if output:
        run_fragments(["save"], namespace)
        namespace["tool_output"] = output
except ValueError as error:
    fail(str(error))
run_report("communication", namespace)
