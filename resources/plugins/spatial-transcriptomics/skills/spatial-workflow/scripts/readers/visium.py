"""Visium and generic AnnData readers.

Visium spatial coordinates come from Space Ranger's spatial/tissue_positions.csv;
``scanpy.read_visium`` reads the whole output directory and places the coordinates
in ``obsm['spatial']``. A bare matrix (.h5 or matrix dir) carries no coordinates,
so the reader reports that rather than inventing a layout.
"""

import os

import numpy as np

from . import finalize


def _platform_meta(adata, platform, **extra):
    meta = {"platform": platform}
    library = adata.uns.get("spatial")
    if isinstance(library, dict) and library:
        meta["library_id"] = list(library.keys())
    meta.update(extra)
    return meta


def read_visium_dir(path, var_names="gene_symbols"):
    import scanpy as sc

    directory = path
    if not os.path.isdir(os.path.join(path, "spatial")) and os.path.isdir(
        os.path.join(path, "filtered_feature_bc_matrix")
    ):
        directory = path
    adata = sc.read_visium(directory, var_names=var_names)
    adata.obs["in_tissue"] = adata.obs.get("in_tissue", 1)
    return finalize(adata, _platform_meta(adata, "visium", input_kind="visium-dir"))


def read_visium_h5(path, var_names="gene_symbols"):
    import h5py
    import scanpy as sc

    adata = sc.read_10x_h5(path, var_names=var_names)
    adata.var_names_make_unique()
    coords = _visium_h5_coords(path, adata.obs_names)
    if coords is None:
        raise ValueError(
            "this .h5 has a matrix but no tissue positions; pass the Space Ranger "
            "output directory (with spatial/tissue_positions.csv) so spots keep their coordinates"
        )
    adata.obsm["spatial"] = coords
    return finalize(adata, _platform_meta(adata, "visium", input_kind="visium-h5"))


def _visium_h5_coords(path, obs_names):
    import h5py

    with h5py.File(path, "r") as handle:
        if "spatial" not in handle:
            return None
        group = handle["spatial"]
        for key in ("tissue_positions", "tissue_positions_list", "positions"):
            if key in group:
                table = np.asarray(group[key])
                if table.ndim == 2 and table.shape[1] >= 5:
                    return table[:, -2:].astype(float)
    return None


def read_10x_h5(path, var_names="gene_symbols"):
    import scanpy as sc

    adata = sc.read_10x_h5(path, var_names=var_names)
    adata.var_names_make_unique()
    return finalize(adata, _platform_meta(adata, "10x-h5", input_kind="10x-h5"))


def read_10x_dir(path, var_names="gene_symbols"):
    import scanpy as sc

    directory = path
    nested = os.path.join(path, "filtered_feature_bc_matrix")
    if os.path.isdir(nested):
        directory = nested
    adata = sc.read_10x_mtx(directory, var_names=var_names)
    adata.var_names_make_unique()
    return finalize(adata, _platform_meta(adata, "10x-matrix-dir", input_kind="10x-matrix-dir"))


def read_anndata(path):
    import scanpy as sc

    adata = sc.read_h5ad(path)
    if "spatial" not in adata.obsm:
        cols = [c for c in ("x", "y") if c in adata.obs]
        if len(cols) == 2:
            adata.obsm["spatial"] = adata.obs[["x", "y"]].to_numpy(dtype=float)
    platform = "anndata"
    if isinstance(adata.uns.get("spatial_platform"), dict):
        platform = adata.uns["spatial_platform"].get("platform", platform)
    return finalize(adata, _platform_meta(adata, platform, input_kind="anndata"))


def read_loom(path):
    import scanpy as sc

    adata = sc.read_loom(path)
    return finalize(adata, _platform_meta(adata, "loom", input_kind="loom"))