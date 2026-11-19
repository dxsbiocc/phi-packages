#!/usr/bin/env python3
"""Thin adapter for sc_py_trajectory: flags -> config globals -> fragments -> JSON report."""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import emit, fail, flag_int, flag_req, flag_str, parse_flags
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
    "group_by": flag_str(flags, "group_by", "leiden"),
    "root_cluster": flag_str(flags, "root_cluster"),
    "n_dcs": flag_int(flags, "n_dcs", 10),
    "seed": flag_int(flags, "seed", 42),
    "plot_tag": "before",
    "emit": emit,
    "tool_input": input_path,
    "tool_output": None,
}

try:
    run_fragments(["setup", "load", "trajectory"], namespace)
    if output:
        run_fragments(["save"], namespace)
        namespace["tool_output"] = output
except ValueError as error:
    fail(str(error))
run_report("trajectory", namespace)