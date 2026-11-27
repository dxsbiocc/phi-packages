#!/usr/bin/env python3
"""Adapter for scatac_tf_deviation: supports both .h5ad (snapatac2/scanpy) and peak x cell TSV/CSV matrices."""

import contextlib
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import emit, fail, flag_int, flag_req, flag_str, info, parse_flags, script_dir
from _config import config_namespace, structural
from _fragments import fragment_path, run_fragments, run_report

flags = parse_flags()
input_path = flag_req(flags, "input")
output = flag_str(flags, "output")

overrides = {
    **structural(flags, input_path, output),
    "meta_path": flag_str(flags, "meta_path"),
    "motif_matches": flag_str(flags, "motif_matches"),
    "peak_gc": flag_str(flags, "peak_gc"),
    "n_background": flag_int(flags, "n_background", 50),
    "tfdev_output": flag_str(flags, "tfdev_output", "regulon/tf_deviations.tsv"),
}
namespace = config_namespace(flags, overrides)

try:
    if input_path.lower().endswith(".h5ad"):
        run_fragments(["setup", "load", "tf-deviation"], namespace)
        if output:
            run_fragments(["save"], namespace)
            namespace["tool_output"] = output
    else:
        namespace.update({
            "np": np,
            "pd": pd,
            "info": info,
            "figure_paths": {},
            "analysis_steps": [],
            "analysis_notes": [],
            "peak_adata": None,
            "adata": None,
            "tfdev_info": None,
        })
        path = fragment_path("tf-deviation")
        with open(path, encoding="utf-8") as fh:
            src = fh.read()
        with contextlib.redirect_stdout(sys.stderr):
            exec(compile(src, path, "exec"), namespace)
except ValueError as error:
    fail(str(error))

run_report("tf_deviation", namespace)
