#!/usr/bin/env python3
"""Thin adapter for sc_py_markers: flags -> config globals -> fragments -> JSON report."""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import emit, fail, flag_bool, flag_int, flag_num, flag_req, flag_str, parse_flags
from _fragments import run_fragments, run_report

flags = parse_flags()
input_path = flag_req(flags, "input")
output = flag_str(flags, "output")

namespace = {
    "data_path": input_path,
    "project_name": flag_str(flags, "project", "phi"),
    "output_object": output,
    "figures_dir": flag_str(flags, "figures_dir"),
    "tables_dir": flag_str(flags, "tables_dir"),
    "min_genes": None,
    "group_by": flag_str(flags, "group_by", "leiden"),
    "test_use": flag_str(flags, "test", "wilcoxon"),
    "logfc_threshold": flag_num(flags, "logfc_threshold", 0.25),
    "min_pct": flag_num(flags, "min_pct", 0.1),
    "only_positive": flag_bool(flags, "only_positive", True),
    "top_n": flag_int(flags, "top_n", 10),
    "plot_heatmap": flag_bool(flags, "heatmap", True),
    "seed": flag_int(flags, "seed", 42),
    "plot_tag": "before",
    "emit": emit,
    "tool_input": input_path,
    "tool_output": None,
}

try:
    run_fragments(["setup", "load", "markers"], namespace)
    if output:
        run_fragments(["save"], namespace)
        namespace["tool_output"] = output
except ValueError as error:
    fail(str(error))
run_report("markers", namespace)