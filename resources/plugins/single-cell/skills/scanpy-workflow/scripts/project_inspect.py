#!/usr/bin/env python3
"""Inspect a single-cell input and report what it contains.

Accepts a .h5ad, a 10x Cell Ranger matrix directory or .h5, a .loom, or a
dense .csv/.tsv count table. Writes an AnnData .h5ad when --output is given.

The scanpy counterpart of seurat-workflow/project_inspect.R; the JSON keys
match so the agent can route on the same report shape.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import (
    detect_kind,
    emit,
    fail,
    flag_int,
    flag_req,
    flag_str,
    parse_flags,
    read_object,
    require_scanpy,
    sorted_levels,
    write_object,
)

flags = parse_flags()
input_path = flag_req(flags, "input")
output = flag_str(flags, "output")
min_cells = flag_int(flags, "min_cells", 0)
min_features = flag_int(flags, "min_features", 0)
var_names = flag_str(flags, "var_names", "gene_symbols")
if var_names not in ("gene_symbols", "gene_ids"):
    fail("--var_names must be gene_symbols or gene_ids")

sc = require_scanpy()
adata = read_object(input_path, var_names=var_names)
if min_cells > 0:
    sc.pp.filter_genes(adata, min_cells=min_cells)
if min_features > 0:
    sc.pp.filter_cells(adata, min_genes=min_features)

kind = detect_kind(input_path)
has_counts = "counts" in adata.layers
value_class = type(adata.X).__name__

interesting = [
    c
    for c in (
        "sample",
        "condition",
        "group",
        "stim",
        "batch",
        "leiden",
        "cell_type",
        "celltype",
        "n_genes_by_counts",
        "total_counts",
        "pct_counts_mt",
    )
    if c in adata.obs.columns
]
metadata = []
for column in interesting:
    values = adata.obs[column]
    if values.dtype.kind in "ifu":
        metadata.append(
            {
                "column": column,
                "kind": "number",
                "min": round(float(values.min()), 4),
                "max": round(float(values.max()), 4),
                "median": round(float(values.median()), 4),
                "missing": int(values.isna().sum()),
            }
        )
    else:
        levels = sorted_levels(values.dropna())
        metadata.append(
            {
                "column": column,
                "kind": "category",
                "distinct": len(levels),
                "levels": levels[:40],
                "missing": int(values.isna().sum()),
            }
        )

ready_for = []
if has_counts:
    ready_for += ["qc", "preprocess"]
if "X_pca" in adata.obsm:
    ready_for.append("integrate")
if "leiden" in adata.obs.columns:
    ready_for.append("markers")
if "rank_genes_groups" in adata.uns:
    ready_for.append("annotate")
if "neighbors" in adata.uns:
    ready_for.append("trajectory")
if has_counts and len(adata.obs.columns) > 0:
    ready_for.append("pseudobulk")

written = None
if output:
    written = write_object(adata, output)

notes = []
if not has_counts:
    notes.append("no counts layer: re-import the raw matrix before QC or pseudobulk")
if value_class in ("csr_matrix", "csc_matrix"):
    notes.append("sparse matrix retained")
if kind == "10x-matrix-dir":
    notes.append("Cell Ranger matrix read with gene symbols from features.tsv")
if "lognorm" in adata.layers:
    notes.append("lognorm snapshot present (pre-scaling expression)")

emit(
    {
        "input": input_path,
        "kind": kind,
        "output": written,
        "shape": {
            "genes": adata.n_vars,
            "cells": adata.n_obs,
            "layers": [k for k in adata.layers.keys() if k is not None],
            "obsm": list(adata.obsm.keys()),
            "uns": list(adata.uns.keys()),
            "snapshot": "lognorm" if "lognorm" in adata.layers else None,
        },
        "counts_present": has_counts,
        "value_class": value_class,
        "metadata": metadata,
        "metadata_columns": list(adata.obs.columns),
        "candidate_steps": ready_for,
        "notes": notes,
    }
)