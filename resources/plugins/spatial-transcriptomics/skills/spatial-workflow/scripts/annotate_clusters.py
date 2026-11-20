#!/usr/bin/env python3
"""Thin adapter for st_py_annotate: flags -> config globals -> fragments -> JSON report.

Without --panel_path and --mapping_path the annotate/assign-labels fragments are
skipped: the notebook falls back to its built-in example panel, but a script run
must not invent labels nobody asked for. A spatial label names the dominant cell
type or region of a spot/bin, not a single pure cell.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import emit, fail, flag_req, flag_str, parse_flags
from _config import config_namespace, structural
from _fragments import run_fragments, run_report

flags = parse_flags()
input_path = flag_req(flags, "input")
output = flag_str(flags, "output")

overrides = structural(flags, input_path, output)
overrides["min_genes"] = None
namespace = config_namespace(flags, overrides)
panel_path = namespace["panel_path"]
mapping_path = namespace["mapping_path"]

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