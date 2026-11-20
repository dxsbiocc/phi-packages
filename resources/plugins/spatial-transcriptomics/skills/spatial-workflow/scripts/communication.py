#!/usr/bin/env python3
"""Thin adapter for st_py_communication: flags -> config globals -> fragments -> JSON report.

Ligand-receptor cell-cell communication (squidpy ligrec) on a clustered,
normalized .h5ad. Runs offline with a curated starter resource unless
--interactions_path supplies a fuller table. liana/CellPhoneDB core are pip-only
and belong in a project environment.
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
    run_fragments(["setup", "load", "spatial-neighbors", "communication"], namespace)
    if output:
        run_fragments(["save"], namespace)
        namespace["tool_output"] = output
except ValueError as error:
    fail(str(error))
run_report("communication", namespace)