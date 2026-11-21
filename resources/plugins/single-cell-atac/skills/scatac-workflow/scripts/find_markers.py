#!/usr/bin/env python3
"""Thin adapter for scatac_markers: flags -> config globals -> fragments -> JSON.

The clustered chain plus per-cluster peaks, marker regions and pairwise differential
accessibility on the peak matrix.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import fail, flag_bool, flag_req, flag_str, parse_flags
from _config import config_namespace, structural
from _fragments import run_fragments, run_report

flags = parse_flags()
input_path = flag_req(flags, "input")
output = flag_str(flags, "output")
save_object = flag_bool(flags, "save_object", True)

namespace = config_namespace(flags, structural(flags, input_path, output))

try:
    run_fragments(
        ["setup", "blacklist", "load", "qc-metrics", "filter", "features", "embed", "cluster", "peaks", "markers"],
        namespace,
    )
    if output and save_object:
        run_fragments(["save"], namespace)
        namespace["tool_output"] = output
except ValueError as error:
    fail(str(error))
run_report("markers", namespace)