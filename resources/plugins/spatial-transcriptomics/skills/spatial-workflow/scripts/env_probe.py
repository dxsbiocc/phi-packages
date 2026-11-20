#!/usr/bin/env python3
"""Report which spatial Python packages this environment provides, and which
spatial workflow modules can therefore run here.

The spatial counterpart of scanpy-workflow/env_probe.py. Phi environments are
read-only and prebuilt, so this script states what the environment already has
and maps each gap to a concrete next step instead of failing mid-analysis. The
model-based spatial tools (cell2location, tangram-sc, stereopy, liana) are
pip-only and cannot enter this pure-conda locked environment; they need a
project environment, never an in-environment install.
"""

import importlib.util
import json
import os
import platform
import sys
from importlib.metadata import PackageNotFoundError, version

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import emit, flag_str, info, parse_flags

flags = parse_flags()
module = flag_str(flags, "module")
write_panel = flag_str(flags, "write_panel")

# (distribution, import name) -> role, modules that rely on it, core, required.
CATALOGUE = [
    ("scanpy", "scanpy", "core single-cell/spatial workflow", [1, 2, 3, 4, 5, 6], True, True),
    ("anndata", "anndata", "AnnData container", [1, 2, 3, 4, 5, 6], True, True),
    ("numpy", "numpy", "numeric core", [1, 2, 3, 4, 5, 6], True, True),
    ("pandas", "pandas", "tables and marker handoff", [1, 2, 3, 4, 5, 6], True, True),
    ("scipy", "scipy", "sparse matrices, statistics, NNLS deconvolution", [1, 2, 4, 5, 6], True, True),
    ("h5py", "h5py", "Stereo-seq GEF / Visium HDF5 reading", [1, 2, 3, 4, 5, 6], True, True),
    ("matplotlib-base", "matplotlib", "figures", [1, 2, 3, 4, 5, 6], True, True),
    ("scikit-learn", "sklearn", "scaling, PCA helpers, neighbours", [2, 3, 4, 5, 6], True, True),
    ("leidenalg", "leidenalg", "Leiden clustering", [2, 3, 4, 5, 6], True, True),
    ("python-igraph", "igraph", "graph backend for Leiden", [2, 3, 4, 5, 6], True, True),
    ("scikit-misc", "skmisc", "loess for seurat_v3 HVG selection", [2, 3, 4, 5, 6], True, True),
    ("squidpy", "squidpy", "spatial graph, SVG, neighbourhood, ligrec", [3, 4, 5, 6], True, True),
    ("spatialdata", "spatialdata", "spatial container (squidpy dependency)", [3, 4, 5, 6], True, True),
    ("spatialdata-io", "spatialdata_io", "spatial format readers", [3, 4, 5, 6], True, False),
    ("omnipath", "omnipath", "ligand-receptor resource (offline curated set used by default)", [5, 6], True, False),
    ("numba", "numba", "JIT acceleration for spatial statistics", [4, 5, 6], True, False),
    ("seaborn", "seaborn", "cross-sample differential heatmap", [6], True, False),
    ("cell2location", "cell2location", "model-based deconvolution (pip-only)", [5, 6], False, False),
    ("tangram-sc", "tangram", "deep-learning deconvolution (pip-only)", [5, 6], False, False),
    ("stereopy", "stereopy", "BGI Stereo-seq toolkit (pip-only)", [5, 6], False, False),
    ("liana", "liana", "CCC consensus framework (pip-only)", [5, 6], False, False),
    ("scvi-tools", "scvi", "deep generative modelling (pip-only)", [], False, False),
]


def probe(dist, import_name, role, modules, core, required):
    installed = importlib.util.find_spec(import_name) is not None
    ver = None
    if installed:
        try:
            ver = version(dist)
        except PackageNotFoundError:
            ver = None
    return {
        "package": dist,
        "role": role,
        "installed": installed,
        "version": ver,
        "core": core,
        "required": required,
        "modules": [str(m) for m in modules],
    }


packages = [probe(*entry) for entry in CATALOGUE]
by_name = {entry["package"]: entry for entry in packages}

modules = []
for number in range(1, 7):
    needed = [entry for entry in CATALOGUE if number in entry[3]]
    missing_required = [
        dist for dist, _, _, _, _, required in needed if required and not by_name[dist]["installed"]
    ]
    missing_optional = [
        dist for dist, _, _, _, core, required in needed
        if not core and not required and not by_name[dist]["installed"]
    ]
    modules.append(
        {
            "module": number,
            "packages_checked": len(needed),
            "missing_required": missing_required,
            "missing_optional": missing_optional,
            "core_ready": len(missing_required) == 0,
        }
    )

selected = (
    [entry for entry in modules if str(entry["module"]) == module]
    if module
    else modules
)

written = None
if write_panel:
    panel = {
        "T_cell": ["CD3D", "CD3E", "TRAC", "IL7R", "LTB"],
        "B_cell": ["MS4A1", "CD79A", "CD79B", "CD19"],
        "Plasma": ["MZB1", "XBP1", "JCHAIN", "SDC1"],
        "NK": ["NKG7", "KLRD1", "GNLY", "PRF1"],
        "Myeloid": ["LYZ", "S100A8", "S100A9", "CTSS"],
        "Dendritic": ["FCER1A", "CLEC10A", "CD1C", "ITGAX"],
        "Fibroblast": ["COL3A1", "PDGFRA", "DCN", "LAMA2"],
        "Epithelial": ["EPCAM", "KRT8", "KRT18", "CDH1"],
        "Endothelial": ["PECAM1", "VWF", "CDH5", "CLDN5"],
    }
    os.makedirs(os.path.dirname(os.path.abspath(write_panel)), exist_ok=True)
    with open(write_panel, "w") as handle:
        json.dump(panel, handle, indent=2)
    info("panel " + write_panel)
    written = write_panel

installed_count = sum(1 for entry in packages if entry["installed"])

notes = []
if installed_count < len(packages):
    notes.append(
        "the plugin environment is read-only and prebuilt: a missing package needs "
        "a project environment (env_request), not an in-environment install"
    )
notes.append(
    "deconvolution runs scipy NNLS in-plugin; cell2location/tangram-sc/stereopy are "
    "pip-only (torch or compiled deps) and need a project environment"
)
notes.append(
    "cell-cell communication runs squidpy ligrec with a curated offline ligand-receptor "
    "set; pass interactions_path for fuller coverage or use liana in a project environment"
)
notes.append(
    "raw reads to counts (Space Ranger / SAW) and Stereo-seq cell segmentation (SAW cellCut) "
    "are upstream wrappers, not notebook cells"
)

emit(
    {
        "environment": "Python " + platform.python_version(),
        "packages_total": len(packages),
        "packages_installed": installed_count,
        "packages": packages,
        "modules": selected,
        "core_ready_modules": [
            str(entry["module"]) for entry in selected if entry["core_ready"]
        ],
        "blocked_modules": [
            str(entry["module"]) for entry in selected if not entry["core_ready"]
        ],
        "panel_written": written,
        "notes": notes,
    }
)