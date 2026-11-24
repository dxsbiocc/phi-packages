#!/usr/bin/env python3
"""Thin adapter for st_py_trajectory: flags -> config globals -> fragments -> JSON report.

Spatial trajectory (Scanpy PAGA + diffusion pseudotime mapped onto tissue
coordinates) and optional gene-program spatial scoring (sc.tl.score_genes on
the unscaled lognorm layer). Zero new dependencies.
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

try:
    run_fragments(["setup", "load", "spatial-neighbors", "trajectory"], namespace)
    if output:
        run_fragments(["save"], namespace)
        namespace["tool_output"] = output
except ValueError as error:
    fail(str(error))
run_report("trajectory", namespace)
