"""Format readers for spatial transcriptomics inputs.

Each reader returns ``(adata, meta)`` where adata is an AnnData with raw counts
in ``X`` (and a ``counts`` layer added by the caller) and spatial coordinates in
``obsm['spatial']``; meta records the platform, bin and coordinate provenance so
the report layer and the agent can audit what was loaded.

Supported platforms
  * Visium        — Space Ranger matrix directory (spatial/ + matrix) or .h5
  * Visium HD     — square-bin or segmented-cell bin tables (.h5 / .csv)
  * Stereo-seq    — SAW GEF (HDF5, square bin or cellbin) and GEM (TSV)

The Stereo-seq GEF/GEM layout follows BGI SAW v8 output. Cell segmentation for
cellbin is produced upstream by SAW cellCut; this reader consumes it, never
re-segments. Where a file's internal layout is not documented, the reader says
so rather than guessing field names.
"""

import os

import numpy as np


def sniff_h5(path):
    """Classify an .h5 by its HDF5 groups (Visium matrix vs Stereo-seq GEF vs HD)."""
    import h5py

    with h5py.File(path, "r") as handle:
        keys = set(handle.keys())
    if "cellBin" in keys:
        return "stereoseq-cellbin"
    if "geneExp" in keys:
        return "stereoseq-gef"
    if "matrix" in keys and "spatial" in keys:
        return "visium-h5"
    if "matrix" in keys:
        return "10x-h5"
    if "filtered_feature_bc_matrix" in keys or "raw_feature_bc_matrix" in keys:
        return "10x-h5"
    return "unknown-h5"


def read_spatial(path, kind, platform="auto", bin_size=None, cellbin=False, var_names="gene_symbols"):
    """Dispatch to the format reader for ``kind``; return (adata, meta)."""
    from . import stereoseq, visium, visium_hd

    if kind == "visium-dir":
        return visium.read_visium_dir(path, var_names=var_names)
    if kind == "visium-h5":
        return visium.read_visium_h5(path, var_names=var_names)
    if kind == "10x-h5":
        return visium.read_10x_h5(path, var_names=var_names)
    if kind == "10x-matrix-dir":
        return visium.read_10x_dir(path, var_names=var_names)
    if kind == "stereoseq-gef":
        return stereoseq.read_gef(path, bin_size=bin_size, cellbin=cellbin)
    if kind == "stereoseq-cellbin":
        return stereoseq.read_gef(path, bin_size=bin_size, cellbin=True)
    if kind == "stereoseq-gem":
        return stereoseq.read_gem(path, bin_size=bin_size, cellbin=cellbin)
    if kind == "anndata":
        return visium.read_anndata(path)
    if kind == "loom":
        return visium.read_loom(path)
    if kind == "dense-table":
        return visium_hd.read_table(path, bin_size=bin_size)
    if kind == "unknown-h5":
        raise ValueError(
            "unrecognised .h5 layout (no matrix/geneExp/cellBin group); expected a "
            "Visium matrix, a Stereo-seq GEF, or a Visium HD bin file"
        )
    if kind == "unknown-directory":
        raise ValueError(
            "--input is a directory with no spatial/ folder or 10x matrix files; "
            "expected a Space Ranger output (spatial/ + matrix) or a Cell Ranger matrix folder"
        )
    raise ValueError(
        "unsupported --input; expected a Visium directory/.h5, a Stereo-seq .gef/.gem, "
        "a Visium HD bin table, or an .h5ad with spatial coordinates"
    )


def finalize(adata, meta, coord_key="spatial"):
    """Make var names unique and record coordinate provenance. Counts layer is added by caller."""
    adata.var_names_make_unique()
    if coord_key not in adata.obsm:
        raise ValueError(
            "reader produced no obsm['" + coord_key + "']; spatial statistics need coordinates"
        )
    meta = dict(meta)
    meta["n_units"] = int(adata.n_obs)
    meta["n_genes"] = int(adata.n_vars)
    meta["coord_key"] = coord_key
    return adata, meta


def bin_coordinates(x, y, bin_size):
    """Aggregate DNB coordinates to square-bin origins: bin_x = (x // N) * N."""
    x = np.asarray(x)
    y = np.asarray(y)
    if bin_size and bin_size > 1:
        x = (x // bin_size) * bin_size
        y = (y // bin_size) * bin_size
    return x, y