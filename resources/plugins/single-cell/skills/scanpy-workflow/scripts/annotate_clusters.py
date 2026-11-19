#!/usr/bin/env python3
"""Thin adapter for sc_py_annotate: flags -> config globals -> fragments -> JSON report.

Without --panel and --mapping the annotate/assign-labels fragments are skipped:
the notebook falls back to its built-in example panel, but a script run must not
invent labels nobody asked for.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import emit, fail, flag_num, flag_req, flag_str, parse_flags
from _fragments import run_fragments, run_report

flags = parse_flags()
input_path = flag_req(flags, "input")
output = flag_str(flags, "output")
panel_path = flag_str(flags, "panel")
mapping_path = flag_str(flags, "mapping")

namespace = {
    "data_path": input_path,
    "project_name": flag_str(flags, "project", "phi"),
    "output_object": output,
    "figures_dir": flag_str(flags, "figures_dir"),
    "min_genes": None,
    "group_by": flag_str(flags, "group_by", "leiden"),
    "label_column": flag_str(flags, "label_column", "cell_type"),
    "annotation_method": flag_str(flags, "method", "panel"),
    "min_score": flag_num(flags, "min_score", 0),
    "panel_path": panel_path,
    "mapping_path": mapping_path,
    "seed": 42,
    "plot_tag": "before",
    "emit": emit,
    "tool_input": input_path,
    "tool_output": None,
}

steps = ["setup", "load"]
if panel_path:
    steps.append("annotate")
if panel_path or mapping_path:
    steps.append("assign-labels")
if output:
    steps.append("save")

try:
    run_fragments(steps, namespace)
    if output:
        namespace["tool_output"] = output
except ValueError as error:
    fail(str(error))
run_report("annotate", namespace)