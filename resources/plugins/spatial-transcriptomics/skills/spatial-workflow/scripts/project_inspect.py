#!/usr/bin/env python3
"""Inspect a spatial input and report what it contains.

Accepts a Visium Space Ranger directory or .h5, a Stereo-seq GEF/GEM (square bin
or cellbin), a Visium HD bin table (.csv/.tsv), or an .h5ad with spatial
coordinates. Writes an AnnData .h5ad when --output is given.

The spatial counterpart of scanpy-workflow/project_inspect.py; the JSON keys
match so the agent can route on the same report shape, plus platform/bin/coord
provenance that only spatial inputs carry.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import (
    detect_kind,
    emit,
    fail,
    flag_bool,
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
platform = flag_str(flags, "platform", "auto")
bin_size = flag_int(flags, "bin_size")
cellbin = flag_bool(flags, "cellbin", False)
min_units = flag_int(flags, "min_units", 0)
min_genes = flag_int(flags, "min_genes", 0)
var_names = flag_str(flags, "var_names", "gene_symbols")
if var_names not in ("gene_symbols", "gene_ids"):
    fail("--var_names must be gene_symbols or gene_ids")

sc = require_scanpy()
adata = read_object(input_path, platform=platform, bin_size=bin_size, cellbin=cellbin, var_names=var_names)
if min_genes > 0:
    sc.pp.filter_genes(adata, min_cells=min_genes)
if min_units > 0:
    sc.pp.filter_cells(adata, min_genes=min_units)

kind = detect_kind(input_path)
meta = adata.uns.get("spatial_platform", {})
has_counts = "counts" in adata.layers
has_coords = "spatial" in adata.obsm
value_class = type(adata.X).__name__

interesting = [
    c
    for c in (
        "sample",
        "condition",
        "group",
        "batch",
        "region",
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
if has_coords:
    ready_for.append("spatial")
if "leiden" in adata.obs.columns:
    ready_for += ["markers", "annotate"]
if "rank_genes_groups" in adata.uns:
    ready_for.append("annotate")
if "lognorm" in adata.layers and "leiden" in adata.obs.columns:
    ready_for += ["spatial-statistics", "communication"]
if has_counts:
    ready_for.append("deconvolve")

written = None
if output:
    written = write_object(adata, output)

notes = []
if not has_counts:
    notes.append("no counts layer: re-import the raw matrix before QC or deconvolution")
if not has_coords:
    notes.append("no obsm['spatial']: spatial statistics need coordinates; pass the Space Ranger directory or a GEF/GEM")
if value_class in ("csr_matrix", "csc_matrix"):
    notes.append("sparse matrix retained")
if meta.get("bin_auto"):
    notes.append(f"bin size auto-selected ({meta.get('bin_size')}); choosing a bin is a biological decision — set bin_size explicitly")
if meta.get("cellbin"):
    notes.append("cellbin read; segmentation was done upstream by SAW cellCut, not here")
if "lognorm" in adata.layers:
    notes.append("lognorm snapshot present (pre-scaling expression)")

emit(
    {
        "input": input_path,
        "kind": kind,
        "output": written,
        "platform": meta,
        "shape": {
            "genes": adata.n_vars,
            "units": adata.n_obs,
            "layers": [k for k in adata.layers.keys() if k is not None],
            "obsm": list(adata.obsm.keys()),
            "obsp": list(adata.obsp.keys()),
            "uns": list(adata.uns.keys()),
            "snapshot": "lognorm" if "lognorm" in adata.layers else None,
        },
        "counts_present": has_counts,
        "coordinates_present": has_coords,
        "value_class": value_class,
        "metadata": metadata,
        "metadata_columns": list(adata.obs.columns),
        "candidate_steps": ready_for,
        "notes": notes,
    }
)