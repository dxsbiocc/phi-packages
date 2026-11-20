#!/usr/bin/env python3
"""Thin adapter for st_py_qc: flags -> config globals -> fragments -> JSON report.

The analysis lives in assets/fragments (the same cells the notebook renders);
this file only translates CLI flags and picks the step order. Spatial units are
spots/bins/cells, not droplets, so there is no doublet gate.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import emit, fail, flag_bool, flag_req, flag_str, parse_flags
from _config import config_namespace, structural
from _fragments import run_fragments, run_report

flags = parse_flags()
input_path = flag_req(flags, "input")
output = flag_str(flags, "output")
save_object = flag_bool(flags, "save_object", True)

namespace = config_namespace(flags, structural(flags, input_path, output))

try:
    run_fragments(["setup", "load", "qc-metrics", "qc-plots"], namespace)
    namespace["plot_tag"] = "after"
    run_fragments(["filter", "qc-plots"], namespace)
    if output and save_object:
        run_fragments(["save"], namespace)
        namespace["tool_output"] = output
except ValueError as error:
    fail(str(error))
run_report("qc", namespace)