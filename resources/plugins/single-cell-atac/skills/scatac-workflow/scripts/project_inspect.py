#!/usr/bin/env python3
"""Inspect a scATAC input and report what it contains (scatac_inspect).

Accepts a fragments.tsv.gz, a position-sorted BAM, a saved snapatac2 .h5ad, or a 10x
peak/tile matrix directory/.h5. Writes an .h5ad when --output is given. The JSON keys
match the other scatac tools so the agent can route on the same report shape, plus the
genome/input-kind provenance that scATAC inputs carry.
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import (
    emit,
    fail,
    flag_req,
    flag_str,
    obs_col,
    obs_columns,
    obs_has,
    parse_flags,
    sorted_levels,
)
from _config import config_namespace, structural
from _fragments import run_fragments

flags = parse_flags()
input_path = flag_req(flags, "input")
output = flag_str(flags, "output")

namespace = config_namespace(flags, structural(flags, input_path, output))

try:
    run_fragments(["setup", "load", "qc-metrics"], namespace)
    if output:
        run_fragments(["save"], namespace)
        namespace["tool_output"] = output
except ValueError as error:
    fail(str(error))

adata = namespace["adata"]
candidates = [
    "sample", "condition", "group", "batch", "region", "leiden", "cell_type", "celltype",
    "n_fragment", "frac_dup", "frac_mito", "tsse",
]
present = obs_columns(adata, candidates)
metadata = []
for column in present:
    values = obs_col(adata, column)
    if values is None:
        continue
    if values.dtype.kind in "ifu":
        v = values.astype(float)
        metadata.append(
            {
                "column": column,
                "kind": "number",
                "min": round(float(v.min()), 4),
                "max": round(float(v.max()), 4),
                "median": round(float(np.median(v)), 4),
            }
        )
    else:
        levels = sorted_levels(values.astype(str))
        metadata.append({"column": column, "kind": "category", "distinct": len(levels), "levels": levels[:40]})

ready_for = []
if obs_has(adata, "n_fragment"):
    ready_for += ["qc", "preprocess"]
if "X_spectral" in namespace.get("embeddings", []):
    ready_for += ["peaks", "markers"]
if obs_has(adata, "leiden"):
    ready_for += ["markers", "annotate"]
if obs_has(adata, "n_fragment"):
    ready_for.append("gene-activity")

emit(
    {
        "input": input_path,
        "kind": namespace.get("loaded_kind"),
        "output": namespace.get("tool_output"),
        "genome": namespace.get("genome_info"),
        "shape": {"cells": adata.n_obs, "features": adata.n_vars},
        "metadata": metadata,
        "metadata_columns": present,
        "candidate_steps": sorted(set(ready_for)),
        "notes": namespace.get("analysis_notes", []),
    }
)