"""Stereo-seq readers for BGI SAW output (GEF HDF5 and GEM TSV).

Square-bin GEF (``/geneExp/binN``)
  * ``gene``       — compound (gene: S64, offset: uint, count: uint); gene i owns
                     expression rows [offset[i], offset[i]+count[i]).
  * ``expression`` — compound (x: uint, y: uint, count: uint); x/y are DNB
                     coordinates, aggregated to a bin origin by (x // N) * N.
  * root attrs carry minX/minY/maxX/maxY/resolution. A tissue.gef holds several
    bins (e.g. 1/5/10/20/50/100/150/200); choosing one is a biological decision.

Cellbin GEF (``/cellBin``)
  * Segmentation is produced upstream by SAW cellCut. ``cellExp`` links cells to
    genes/counts; ``gene`` maps gene IDs to names. Internal field names vary by
    SAW version, so the reader inspects the compound dtype and adapts, reporting
    the fields it used rather than assuming a fixed layout.

GEM (SAW v8 TSV)
  * ``#`` comment lines precede a header. Square bin: geneID, geneName, x, y,
    MIDCount, ExonCount. Cellbin adds a CellID column. Tab-separated.
"""

import gzip
import os

import numpy as np
import pandas as pd
from scipy import sparse

from . import bin_coordinates, finalize


def _open_text(path):
    if path.lower().endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def list_gef_bins(path):
    import h5py

    with h5py.File(path, "r") as handle:
        if "geneExp" not in handle:
            return []
        bins = []
        for key in handle["geneExp"].keys():
            if key.startswith("bin"):
                try:
                    bins.append(int(key[3:]))
                except ValueError:
                    continue
    return sorted(bins)


def choose_bin(available, bin_size):
    """Pick a bin deterministically: explicit choice > the only bin > ~10µm default."""
    if not available:
        raise ValueError("GEF has no /geneExp/binN group; expected a square-bin Stereo-seq GEF")
    if bin_size is not None:
        if bin_size not in available:
            raise ValueError(
                f"bin {bin_size} not in this GEF; available bins: {available}"
            )
        return bin_size, False
    if len(available) == 1:
        return available[0], False
    for preferred in (20, 50, 10, 100):
        if preferred in available:
            return preferred, True
    return min(available), True


def read_gef(path, bin_size=None, cellbin=False):
    import h5py

    with h5py.File(path, "r") as handle:
        if cellbin or "cellBin" in handle:
            return _read_cellbin(handle, path)
        return _read_square_bin(handle, path, bin_size)


def _read_square_bin(handle, path, bin_size):
    available = []
    for key in handle["geneExp"].keys():
        if key.startswith("bin"):
            try:
                available.append(int(key[3:]))
            except ValueError:
                continue
    chosen, auto = choose_bin(sorted(available), bin_size)
    group = handle["geneExp"][f"bin{chosen}"]
    gene = np.asarray(group["gene"])
    expression = np.asarray(group["expression"])

    genes = np.array([g.decode() if isinstance(g, bytes) else str(g) for g in gene["gene"]])
    offsets = gene["offset"].astype(np.int64)
    counts = gene["count"].astype(np.int64)

    x = expression["x"].astype(np.int64)
    y = expression["y"].astype(np.int64)
    mid = expression["count"].astype(np.float64)
    bx, by = bin_coordinates(x, y, chosen)

    bins, bin_index = np.unique(np.stack([bx, by], axis=1), axis=0, return_inverse=True)
    n_bins = bins.shape[0]
    n_genes = genes.shape[0]

    rows = bin_index
    cols = np.repeat(np.arange(n_genes, dtype=np.int64), counts)
    data = mid
    matrix = sparse.coo_matrix((data, (rows, cols)), shape=(n_bins, n_genes)).tocsr()
    matrix.sum_duplicates()

    import anndata as ad

    obs = pd.DataFrame(index=[f"{int(a)}_{int(b)}" for a, b in bins])
    var = pd.DataFrame(index=pd.Index(genes, name=None))
    adata = ad.AnnData(X=matrix, obs=obs, var=var)
    adata.obsm["spatial"] = bins.astype(float)
    attrs = dict(handle.attrs)
    meta = {
        "platform": "stereoseq",
        "input_kind": "stereoseq-gef",
        "bin_size": chosen,
        "bin_auto": auto,
        "available_bins": sorted(available),
        "resolution": attrs.get("resolution"),
        "bin_origin": "DNB (x//N)*N",
    }
    return finalize(adata, meta)


def _read_cellbin(handle, path):
    group = handle["cellBin"]
    fields = {name: group[name].dtype.names for name in group.keys() if hasattr(group[name], "dtype")}
    if "cellExp" not in group or "gene" not in group:
        raise ValueError(
            "cellbin GEF missing cellExp/gene groups; found " + ", ".join(group.keys())
        )
    gene = np.asarray(group["gene"])
    gene_fields = gene.dtype.names
    name_field = next((f for f in ("gene", "geneName", "name") if f in gene_fields), gene_fields[0])
    genes = np.array([g.decode() if isinstance(g, bytes) else str(g) for g in gene[name_field]])
    gene_id_field = next((f for f in ("geneID", "id", "offset") if f in gene_fields), None)

    cell_exp = np.asarray(group["cellExp"])
    exp_fields = cell_exp.dtype.names
    cell_field = next((f for f in ("cellID", "cell", "id") if f in exp_fields), None)
    count_field = next((f for f in ("MIDCount", "count", "umi", "exp") if f in exp_fields), None)
    gene_ref_field = next((f for f in ("geneID", "gene", "geneIndex") if f in exp_fields), None)
    if cell_field is None or count_field is None or gene_ref_field is None:
        raise ValueError(
            "cellbin cellExp fields not recognised: " + ", ".join(exp_fields or ())
            + "; expected cell/gene/count columns. SAW cellbin layouts vary by version."
        )

    cell_ids = cell_exp[cell_field].astype(np.int64)
    gene_ref = cell_exp[gene_ref_field].astype(np.int64)
    counts = cell_exp[count_field].astype(np.float64)

    cells, cell_index = np.unique(cell_ids, return_inverse=True)
    if gene_id_field is not None:
        gene_ids = gene[gene_id_field].astype(np.int64)
        order = np.argsort(gene_ids)
        gene_pos = np.searchsorted(gene_ids[order], gene_ref)
        cols = order[np.clip(gene_pos, 0, len(genes) - 1)]
    else:
        cols = np.clip(gene_ref, 0, len(genes) - 1)

    matrix = sparse.coo_matrix(
        (counts, (cell_index, cols)), shape=(cells.shape[0], genes.shape[0])
    ).tocsr()
    matrix.sum_duplicates()

    coords = _cellbin_coords(group, cells)
    import anndata as ad

    adata = ad.AnnData(
        X=matrix,
        obs=pd.DataFrame(index=[str(int(c)) for c in cells]),
        var=pd.DataFrame(index=pd.Index(genes)),
    )
    if coords is not None:
        adata.obsm["spatial"] = coords
    else:
        raise ValueError(
            "cellbin GEF has no cell centroid coordinates (cell/cellBorder layout not recognised); "
            "spatial statistics need centroids — fields seen: " + ", ".join(group.keys())
        )
    meta = {
        "platform": "stereoseq",
        "input_kind": "stereoseq-cellbin",
        "bin_size": None,
        "cellbin": True,
        "segmentation": "upstream SAW cellCut",
        "cellExp_fields": list(exp_fields or ()),
    }
    return finalize(adata, meta)


def _cellbin_coords(group, cells):
    for name in ("cell", "cellBorder", "cellCentroid"):
        if name not in group:
            continue
        table = np.asarray(group[name])
        if table.dtype.names is None:
            continue
        fields = table.dtype.names
        id_field = next((f for f in ("cellID", "cell", "id") if f in fields), None)
        x_field = next((f for f in ("x", "centerX", "cx") if f in fields), None)
        y_field = next((f for f in ("y", "centerY", "cy") if f in fields), None)
        if id_field and x_field and y_field:
            frame = pd.DataFrame({
                "id": table[id_field].astype(np.int64),
                "x": table[x_field].astype(float),
                "y": table[y_field].astype(float),
            }).drop_duplicates("id").set_index("id")
            aligned = frame.reindex(cells)
            if aligned[["x", "y"]].notna().all().all():
                return aligned[["x", "y"]].to_numpy()
    return None


def read_gem(path, bin_size=None, cellbin=False):
    frame = _read_gem_frame(path)
    if cellbin or "CellID" in frame.columns:
        return _gem_cellbin(frame)
    return _gem_square_bin(frame, bin_size)


def _read_gem_frame(path):
    with _open_text(path) as handle:
        header = None
        for line in handle:
            if line.startswith("#"):
                continue
            header = line.rstrip("\n").split("\t")
            break
        if header is None:
            raise ValueError("GEM file has no header line")
        rows = []
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            rows.append(line.rstrip("\n").split("\t"))
    frame = pd.DataFrame(rows, columns=header)
    return frame


def _gem_square_bin(frame, bin_size):
    required = {"geneID", "geneName", "x", "y", "MIDCount"}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError("GEM missing columns: " + ", ".join(sorted(missing)))
    genes = frame["geneName"].astype(str).to_numpy()
    x = frame["x"].astype(np.int64).to_numpy()
    y = frame["y"].astype(np.int64).to_numpy()
    mid = frame["MIDCount"].astype(np.float64).to_numpy()

    chosen = bin_size
    auto = False
    if chosen is None:
        chosen = 1
    bx, by = bin_coordinates(x, y, chosen)
    bins, bin_index = np.unique(np.stack([bx, by], axis=1), axis=0, return_inverse=True)
    gene_names, gene_index = np.unique(genes, return_inverse=True)

    matrix = sparse.coo_matrix(
        (mid, (bin_index, gene_index)), shape=(bins.shape[0], gene_names.shape[0])
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
        "platform": "stereoseq",
        "input_kind": "stereoseq-gem",
        "bin_size": chosen,
        "bin_auto": auto,
        "bin_origin": "DNB (x//N)*N",
    }
    return finalize(adata, meta)


def _gem_cellbin(frame):
    required = {"CellID", "geneName", "x", "y", "MIDCount"}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError("cellbin GEM missing columns: " + ", ".join(sorted(missing)))
    cells, cell_index = np.unique(frame["CellID"].astype(np.int64).to_numpy(), return_inverse=True)
    gene_names, gene_index = np.unique(frame["geneName"].astype(str).to_numpy(), return_inverse=True)
    mid = frame["MIDCount"].astype(np.float64).to_numpy()
    matrix = sparse.coo_matrix(
        (mid, (cell_index, gene_index)), shape=(cells.shape[0], gene_names.shape[0])
    ).tocsr()
    matrix.sum_duplicates()

    coords = frame.groupby(frame["CellID"].astype(np.int64))[["x", "y"]].first()
    coords = coords.reindex(cells).to_numpy(dtype=float)

    import anndata as ad

    adata = ad.AnnData(
        X=matrix,
        obs=pd.DataFrame(index=[str(int(c)) for c in cells]),
        var=pd.DataFrame(index=pd.Index(gene_names)),
    )
    adata.obsm["spatial"] = coords
    meta = {
        "platform": "stereoseq",
        "input_kind": "stereoseq-gem-cellbin",
        "bin_size": None,
        "cellbin": True,
        "segmentation": "upstream SAW cellCut",
    }
    return finalize(adata, meta)