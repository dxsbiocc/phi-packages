#!/usr/bin/env python3
"""Adapter for scatac_peak2gene: supports both .h5ad and peak/gene TSV/CSV matrices."""

import contextlib
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import emit, fail, flag_int, flag_req, flag_str, info, parse_flags
from _config import config_namespace, structural
from _fragments import fragment_path, run_fragments, run_report

flags = parse_flags()
input_path = flag_req(flags, "input")

overrides = {
    **structural(flags, input_path, None),
    "meta_path": flag_str(flags, "meta_path"),
    "rna_input": flag_req(flags, "rna_input"),
    "rna_meta": flag_str(flags, "rna_meta"),
    "gene_tss": flag_str(flags, "gene_tss"),
    "peak2gene_max_dist": flag_int(flags, "max_distance", 250000),
    "peak2gene_n_null": flag_int(flags, "n_null", 200),
    "peak2gene_output": flag_str(flags, "links_output", "peak2gene/peak2gene_links.tsv"),
}
namespace = config_namespace(flags, overrides)

try:
    if input_path.lower().endswith(".h5ad"):
        run_fragments(["setup", "load", "peak2gene"], namespace)
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
            "peak2gene_info": None,
        })
        path = fragment_path("peak2gene")
        with open(path, encoding="utf-8") as fh:
            src = fh.read()
        with contextlib.redirect_stdout(sys.stderr):
            exec(compile(src, path, "exec"), namespace)
except ValueError as error:
    fail(str(error))

run_report("peak2gene", namespace)
