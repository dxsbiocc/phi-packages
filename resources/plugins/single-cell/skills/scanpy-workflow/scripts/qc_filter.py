#!/usr/bin/env python3
"""Thin adapter for sc_py_qc: flags -> config globals -> fragments -> JSON report.

The analysis lives in assets/fragments (the same cells the notebook
renders); this file only translates CLI flags and picks the step order.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import emit, fail, flag_bool, flag_int, flag_num, flag_nums, flag_req, flag_str, parse_flags
from _fragments import run_fragments, run_report

flags = parse_flags()
input_path = flag_req(flags, "input")
output = flag_str(flags, "output")
save_object = flag_bool(flags, "save_object", True)

namespace = {
    "data_path": input_path,
    "project_name": flag_str(flags, "project", "phi"),
    "output_object": output,
    "figures_dir": flag_str(flags, "figures_dir"),
    "min_genes": flag_num(flags, "feature_low"),
    "max_genes": flag_num(flags, "feature_high"),
    "min_counts": flag_num(flags, "count_low"),
    "max_counts": flag_num(flags, "count_high"),
    "mt_max": flag_num(flags, "mt_max", 10),
    "mt_pattern": flag_str(flags, "mt_pattern", "^MT-"),
    "ribo_pattern": flag_str(flags, "ribo_pattern", "^RP[SL]"),
    "hb_pattern": flag_str(flags, "hb_pattern", "^HB[ABDEGMQZ][0-9]*$"),
    "quantiles": flag_nums(flags, "quantiles", [0.01, 0.99]),
    "sample_column": flag_str(flags, "sample_column", "sample"),
    "detect_doublets": flag_bool(flags, "detect_doublets", False),
    "doublet_expected_rate": flag_num(flags, "doublet_expected_rate"),
    "seed": flag_int(flags, "seed", 42),
    "plot_tag": "before",
    "emit": emit,
    "tool_input": input_path,
    "tool_output": None,
}

try:
    run_fragments(["setup", "load", "qc-metrics", "doublets", "qc-plots"], namespace)
    namespace["plot_tag"] = "after"
    run_fragments(["filter", "qc-plots"], namespace)
    if output and save_object:
        run_fragments(["save"], namespace)
        namespace["tool_output"] = output
except ValueError as error:
    fail(str(error))
run_report("qc", namespace)