"""Visium HD and generic spatial bin tables.

Visium HD square bins and segmented cells are exported as tables. The common
long format has one row per (bin, gene) with coordinate and count columns; a wide
format has genes x bins. Column names vary by export, so the reader matches
common aliases and reports what it used. A table without coordinates cannot
support spatial statistics and is rejected with guidance.
"""

import numpy as np
import pandas as pd
from scipy import sparse

from . import finalize

_GENE = ("gene", "gene_name", "genename", "geneid", "gene_id", "feature", "symbol")
_X = ("x", "x_coord", "xcoord", "px", "bin_x", "center_x", "centroid_x")
_Y = ("y", "y_coord", "ycoord", "py", "bin_y", "center_y", "centroid_y")
_COUNT = ("count", "counts", "midcount", "mid_count", "umi", "umis", "value", "expression")
_BIN = ("barcode", "bin", "bin_id", "binid", "cell_id", "cellid", "spot", "unit")


def _find(columns, aliases):
    lower = {str(c).lower(): c for c in columns}
    for alias in aliases:
        if alias in lower:
            return lower[alias]
    return None


def _read_frame(path):
    lower = path.lower()
    sep = "," if lower.endswith((".csv", ".csv.gz")) else "\t"
    return pd.read_csv(path, sep=sep)


def read_table(path, bin_size=None):
    frame = _read_frame(path)
    gene_col = _find(frame.columns, _GENE)
    x_col = _find(frame.columns, _X)
    y_col = _find(frame.columns, _Y)
    count_col = _find(frame.columns, _COUNT)

    if gene_col and x_col and y_col and count_col:
        return _long_format(frame, gene_col, x_col, y_col, count_col, bin_size)
    return _wide_format(frame, path)


def _long_format(frame, gene_col, x_col, y_col, count_col, bin_size):
    genes = frame[gene_col].astype(str).to_numpy()
    x = frame[x_col].astype(np.int64).to_numpy()
    y = frame[y_col].astype(np.int64).to_numpy()
    counts = frame[count_col].astype(np.float64).to_numpy()
    if bin_size and bin_size > 1:
        x = (x // bin_size) * bin_size
        y = (y // bin_size) * bin_size
    bins, bin_index = np.unique(np.stack([x, y], axis=1), axis=0, return_inverse=True)
    gene_names, gene_index = np.unique(genes, return_inverse=True)
    matrix = sparse.coo_matrix(
        (counts, (bin_index, gene_index)), shape=(bins.shape[0], gene_names.shape[0])
    ).tocsr()
    matrix.sum_duplicates()

    import anndata as ad

    adata = ad.AnnData(
        X=matrix,
        obs=pd.DataFrame(index=[f"{int(a)}_{int(b)}" for a, b in bins]),
        var=pd.DataFrame(index=pd.Index(gene_names)),
    )
    adata.obsm["spatial"] = bins.astype(float)
    meta = {
        "platform": "visium-hd",
        "input_kind": "spatial-table-long",
        "bin_size": bin_size,
        "columns": {"gene": gene_col, "x": x_col, "y": y_col, "count": count_col},
    }
    return finalize(adata, meta)


def _wide_format(frame, path):
    x_col = _find(frame.columns, _X)
    y_col = _find(frame.columns, _Y)
    bin_col = _find(frame.columns, _BIN)
    gene_col = _find(frame.columns, _GENE)

    if gene_col is not None:
        frame = frame.set_index(gene_col)
        x_col = _find(frame.columns, _X)
        y_col = _find(frame.columns, _Y)
        if x_col and y_col:
            coords = frame[[x_col, y_col]].to_numpy(dtype=float)
            matrix = frame.drop(columns=[x_col, y_col]).to_numpy(dtype=float)
            import anndata as ad

            adata = ad.AnnData(
                X=sparse.csr_matrix(matrix.T),
                obs=pd.DataFrame(index=[str(c) for c in frame.drop(columns=[x_col, y_col]).columns]),
                var=pd.DataFrame(index=pd.Index(frame.index.astype(str))),
            )
            adata.obsm["spatial"] = coords
            return finalize(adata, {"platform": "visium-hd", "input_kind": "spatial-table-wide"})

    raise ValueError(
        "table has no recognised gene/x/y/count columns; expected a long spatial bin "
        "table (gene, x, y, count) or a wide genes-x-bins matrix with coordinate columns. "
        "Columns seen: " + ", ".join(str(c) for c in frame.columns)
    )