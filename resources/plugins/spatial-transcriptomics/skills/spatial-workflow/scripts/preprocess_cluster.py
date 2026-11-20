#!/usr/bin/env python3
"""Thin adapter for st_py_preprocess: flags -> config globals -> fragments -> JSON report.

Normalize, variable genes, scaled PCA, Leiden clusters and a UMAP in expression
space, then the spatial adjacency graph and cluster maps in tissue space. Reads
an already-QC'd .h5ad (min_genes gate disabled so filtering is not repeated).
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
    run_fragments(
        ["setup", "load", "normalize", "hvg", "scale", "pca", "cluster",
         "spatial-neighbors", "spatial-plots"],
        namespace,
    )
    if output:
        run_fragments(["save"], namespace)
        namespace["tool_output"] = output
except ValueError as error:
    fail(str(error))
run_report("preprocess", namespace)