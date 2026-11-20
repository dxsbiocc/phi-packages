#!/usr/bin/env python3
"""Report which single-cell Python packages this environment provides, and
which NBIS workshop labs can therefore run here.

The scanpy counterpart of seurat-workflow/env_probe.R. The NBIS labs install
their own packages on every run; Phi environments are read-only and prebuilt,
so this script states what the environment already has and maps each gap to
a concrete next step instead of failing mid-analysis.
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

# (distribution, import name) -> role, and which NBIS labs rely on it.
CATALOGUE = [
    ("scanpy", "scanpy", "core single-cell workflow", [1, 2, 3, 4, 5, 6, 7], True, False),
    ("anndata", "anndata", "AnnData container", [1, 2, 3, 4, 5, 6, 7], True, False),
    ("numpy", "numpy", "numeric core", [1, 2, 3, 4, 5, 6, 7], True, False),
    ("pandas", "pandas", "tables and pseudobulk handoff", [1, 5], True, False),
    ("scipy", "scipy", "sparse matrices and statistics", [1, 2, 4], True, False),
    ("matplotlib-base", "matplotlib", "figures", [1, 2, 3, 4, 6, 7], True, False),
    ("scikit-learn", "sklearn", "scaling, PCA helpers, neighbours", [2], True, False),
    ("leidenalg", "leidenalg", "Leiden clustering", [2, 3], True, False),
    ("python-igraph", "igraph", "graph backend for Leiden and PAGA", [2, 3, 7], True, False),
    ("scikit-misc", "skmisc", "loess for seurat_v3 HVG selection", [2], True, False),
    ("scrublet", "scrublet", "doublet detection (lab 1)", [1], False, False),
    ("harmonypy", "harmonypy", "batch correction (lab 3)", [3], False, True),
    ("bbknn", "bbknn", "batch-balanced KNN (lab 3 alternative)", [3], False, False),
    ("scanorama", "scanorama", "panorama integration (lab 3 alternative)", [3], False, False),
    ("celltypist", "celltypist", "reference-based annotation (lab 6)", [6], False, False),
    ("scvi-tools", "scvi", "deep generative integration and DE", [3, 4], False, False),
    ("gseapy", "gseapy", "pathway enrichment on marker tables", [4], False, False),
    ("squidpy", "squidpy", "spatial transcriptomics", [], False, False),
]


def probe(dist, import_name, role, labs, core, required):
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
        "modules": [str(lab) for lab in labs],
    }


packages = [probe(*entry) for entry in CATALOGUE]
by_name = {entry["package"]: entry for entry in packages}

modules = []
for number in range(1, 8):
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
        "T_cell": ["CD3D", "CD3E", "TRAC", "TRBC1", "IL7R", "LTB"],
        "B_cell": ["MS4A1", "CD79A", "CD79B", "CD74", "HLA-DRA", "CD19"],
        "Plasma": ["MZB1", "XBP1", "JCHAIN", "SDC1", "TNFRSF17", "IGKC"],
        "NK": ["NKG7", "KLRD1", "GNLY", "PRF1", "TRDC"],
        "Myeloid": ["LYZ", "S100A8", "S100A9", "LGALS3", "CTSS", "MS4A7"],
        "Dendritic": ["FCER1A", "CLEC10A", "CD1C", "ITGAX", "LILRA4", "IRF7"],
        "Fibroblast": ["COL3A1", "PDGFRA", "LAMA2", "FBLN2", "DCN"],
        "Epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1"],
        "Endothelial": ["PECAM1", "VWF", "CDH5", "CLDN5"],
        "Erythrocyte": ["HBB", "HBA1", "HBA2", "ALAS2", "SLC4A1", "GYPA"],
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
if not by_name["bbknn"]["installed"]:
    notes.append(
        "bbknn has no python 3.13 build on bioconda; integration runs harmony or "
        "combat instead (sc_py_integrate)"
    )
notes.append("the NBIS labs install packages on every run; these scripts never do")

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