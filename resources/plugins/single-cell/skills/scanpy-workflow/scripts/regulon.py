#!/usr/bin/env python3
"""Thin adapter for sc_py_regulon: flags -> config globals -> fragments -> JSON report."""

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
    "group_by": flag_str(flags, "group_by", "leiden"),
    "regulons_path": flag_str(flags, "regulons_path"),
    "regulon_method": flag_str(flags, "regulon_method", "aucell"),
    "regulon_top_frac": flag_num(flags, "top_frac", 0.05),
    "regulon_min_targets": flag_int(flags, "min_targets", 5),
    "regulon_output": flag_str(flags, "regulon_output", "regulon/regulon_activity.tsv"),
    "seed": flag_int(flags, "seed", 42),
    "plot_tag": "before",
    "emit": emit,
    "tool_input": input_path,
    "tool_output": None,
}

try:
    run_fragments(["setup", "load", "regulon"], namespace)
    if output:
        run_fragments(["save"], namespace)
        namespace["tool_output"] = output
except ValueError as error:
    fail(str(error))
run_report("regulon", namespace)
