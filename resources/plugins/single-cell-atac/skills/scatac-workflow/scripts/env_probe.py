#!/usr/bin/env python3
"""Report which scATAC Python packages this environment provides, and which
scATAC workflow modules can therefore run here.

Phi environments are read-only and prebuilt, so this script states what the
environment already has and maps each gap to a concrete next step instead of
failing mid-analysis. snapatac2 is the single engine and covers the whole
in-plugin chain (fragment import, QC/TSSe, scrublet, spectral embedding,
Leiden/UMAP, MACS3 peaks, marker regions, differential accessibility, gene
activity, motif enrichment, harmony/scanorama integration). The R-route suites
(ArchR, Signac, chromVAR) and pip-only tools (episcanpy under numpy 2, scvi-tools,
magic-impute) cannot enter this pure-conda locked environment; they need a
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
# macs3 ships its Python module as MACS3 (uppercase); snapatac2.macs3 drives it.
CATALOGUE = [
    ("snapatac2", "snapatac2", "core scATAC engine (Rust file-backed AnnData)", [1, 2, 3, 4, 5, 6], True, True),
    ("scanpy", "scanpy", "gene-activity downstream, panel scoring, scRNA reference exchange", [4, 5, 6], True, True),
    ("anndata", "anndata", "AnnData container", [1, 2, 3, 4, 5, 6], True, True),
    ("numpy", "numpy", "numeric core", [1, 2, 3, 4, 5, 6], True, True),
    ("pandas", "pandas", "tables and marker handoff", [1, 2, 3, 4, 5, 6], True, True),
    ("polars", "polars", "snapatac2 columnar backend, peak/diff tables", [1, 2, 3, 4, 5, 6], True, True),
    ("pyarrow", "pyarrow", "fragment / peak matrix IO", [1, 2, 3, 4, 5, 6], True, True),
    ("scipy", "scipy", "sparse matrices and statistics", [2, 3, 4, 5, 6], True, True),
    ("h5py", "h5py", "h5ad reading and writing", [1, 2, 3, 4, 5, 6], True, True),
    ("scikit-learn", "sklearn", "spectral helpers, neighbours, reference kNN transfer", [2, 3, 4, 5, 6], True, True),
    ("leidenalg", "leidenalg", "Leiden clustering", [2, 3, 4, 5, 6], True, True),
    ("python-igraph", "igraph", "graph backend for Leiden", [2, 3, 4, 5, 6], True, True),
    ("umap-learn", "umap", "UMAP embedding", [2, 3, 4, 5, 6], True, True),
    ("macs3", "MACS3", "per-cluster peak calling", [3, 4, 5, 6], True, True),
    ("plotly", "plotly", "snapatac2 figures (HTML/SVG native)", [1, 2, 3, 4, 5, 6], True, True),
    ("matplotlib-base", "matplotlib", "scanpy gene-activity / annotation figures", [4, 5, 6], True, True),
    ("kaleido", "kaleido", "static PNG/PDF export of plotly figures (HTML works without it)", [1, 2, 3, 4, 5, 6], True, False),
    ("episcanpy", "episcanpy", "alternative scATAC toolkit (bioconda, but numpy 2 incompatible)", [], False, False),
    ("chromvar", "chromvar", "per-cell TF deviation scoring (R/Bioconductor)", [], False, False),
    ("ArchR", "ArchR", "R scATAC suite (R/Bioconductor)", [], False, False),
    ("Signac", "Signac", "R/Seurat scATAC extension (R/Bioconductor)", [], False, False),
    ("scvi-tools", "scvi", "deep generative modelling (pip-only, torch)", [], False, False),
    ("magic-impute", "magic", "diffusion imputation (pip-only)", [], False, False),
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
    "snapatac2 is the single engine and covers the whole in-plugin chain; the TF layer "
    "is motif_enrichment on marker regions, not chromVAR per-cell deviations"
)
notes.append(
    "chromVAR / ArchR / Signac are R-route suites and episcanpy breaks under numpy 2; "
    "scvi-tools / magic-impute are pip-only — all need a project environment, never an "
    "in-environment install"
)
notes.append(
    "motif enrichment needs a genome FASTA (genome_fasta) and gene activity / TSSe need a "
    "gene annotation (gene_anno); hg38/mm10/hg19/mm39 chrom_sizes are bundled, large per-genome "
    "FASTA/GFF3 are user input-path or a reported download gap"
)
notes.append(
    "raw reads to fragments (cellranger-atac / chromap) and bulk ChIP/CUT&Tag peak calling "
    "are upstream wrappers or the sibling chromatin-profiling plugin, not notebook cells"
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