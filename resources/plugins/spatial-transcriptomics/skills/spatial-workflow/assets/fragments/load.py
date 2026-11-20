# data_path, platform, bin_size and cellbin are set in the config cell above.
# Spatial coordinates land in adata.obsm['spatial']; every later spatial statistic
# reads them from there. Raw counts are kept in adata.layers['counts'].
import gzip

import h5py
from scipy import sparse


def _bin_xy(x, y, n):
    x = np.asarray(x, dtype=np.int64)
    y = np.asarray(y, dtype=np.int64)
    if n and n > 1:
        x = (x // n) * n
        y = (y // n) * n
    return x, y


def _bins_to_adata(bx, by, genes, gene_index, counts):
    bins, bin_index = np.unique(np.stack([bx, by], axis=1), axis=0, return_inverse=True)
    gene_names, gene_cols = np.unique(genes, return_inverse=True)
    matrix = sparse.coo_matrix(
        (counts, (bin_index, gene_cols)), shape=(bins.shape[0], gene_names.shape[0])
    ).tocsr()
    matrix.sum_duplicates()
    out = sc.AnnData(
        X=matrix,
        obs=pd.DataFrame(index=[f"{int(a)}_{int(b)}" for a, b in bins]),
        var=pd.DataFrame(index=pd.Index(gene_names)),
    )
    out.obsm["spatial"] = bins.astype(float)
    return out


def _read_gef_square(path, chosen):
    with h5py.File(path, "r") as handle:
        available = sorted(int(k[3:]) for k in handle["geneExp"] if k.startswith("bin"))
        auto = False
        if chosen is None:
            if len(available) == 1:
                chosen = available[0]
            else:
                chosen = next((b for b in (20, 50, 10, 100) if b in available), min(available))
                auto = True
        if chosen not in available:
            raise ValueError(f"bin {chosen} not in GEF; available: {available}")
        group = handle["geneExp"][f"bin{chosen}"]
        gene = np.asarray(group["gene"])
        expression = np.asarray(group["expression"])
        genes = np.array([g.decode() if isinstance(g, bytes) else str(g) for g in gene["gene"]])
        offsets = gene["offset"].astype(np.int64)
        per_gene = gene["count"].astype(np.int64)
        bx, by = _bin_xy(expression["x"], expression["y"], chosen)
        bins, bin_index = np.unique(np.stack([bx, by], axis=1), axis=0, return_inverse=True)
        rows = bin_index
        cols = np.repeat(np.arange(genes.shape[0], dtype=np.int64), per_gene)
        matrix = sparse.coo_matrix(
            (expression["count"].astype(np.float64), (rows, cols)),
            shape=(bins.shape[0], genes.shape[0]),
        ).tocsr()
        matrix.sum_duplicates()
        out = sc.AnnData(
            X=matrix,
            obs=pd.DataFrame(index=[f"{int(a)}_{int(b)}" for a, b in bins]),
            var=pd.DataFrame(index=pd.Index(genes)),
        )
        out.obsm["spatial"] = bins.astype(float)
        info(f"Stereo-seq GEF square bin {chosen} (auto={auto}); bins available: {available}")
        return out, {"platform": "stereoseq", "bin_size": chosen, "bin_auto": auto,
                     "available_bins": available, "input_kind": "stereoseq-gef"}


def _read_gem(path, chosen, as_cellbin):
    opener = gzip.open if path.lower().endswith(".gz") else open
    with opener(path, "rt") as handle:
        header = None
        for line in handle:
            if not line.startswith("#"):
                header = line.rstrip("\n").split("\t")
                break
        rows = [ln.rstrip("\n").split("\t") for ln in handle if ln.strip() and not ln.startswith("#")]
    frame = pd.DataFrame(rows, columns=header)
    if as_cellbin or "CellID" in frame.columns:
        cells, cell_index = np.unique(frame["CellID"].astype(np.int64), return_inverse=True)
        gene_names, gene_index = np.unique(frame["geneName"].astype(str), return_inverse=True)
        matrix = sparse.coo_matrix(
            (frame["MIDCount"].astype(np.float64), (cell_index, gene_index)),
            shape=(cells.shape[0], gene_names.shape[0]),
        ).tocsr()
        matrix.sum_duplicates()
        coords = frame.groupby(frame["CellID"].astype(np.int64))[["x", "y"]].first().reindex(cells)
        out = sc.AnnData(X=matrix, obs=pd.DataFrame(index=[str(int(c)) for c in cells]),
                         var=pd.DataFrame(index=pd.Index(gene_names)))
        out.obsm["spatial"] = coords.to_numpy(dtype=float)
        return out, {"platform": "stereoseq", "cellbin": True, "bin_size": None,
                     "input_kind": "stereoseq-gem-cellbin", "segmentation": "upstream SAW cellCut"}
    chosen = chosen or 1
    bx, by = _bin_xy(frame["x"], frame["y"], chosen)
    out = _bins_to_adata(bx, by, frame["geneName"].astype(str).to_numpy(),
                         None, frame["MIDCount"].astype(np.float64).to_numpy())
    return out, {"platform": "stereoseq", "bin_size": chosen, "input_kind": "stereoseq-gem"}


lower = data_path.lower()
if os.path.isdir(data_path):
    if os.path.isdir(os.path.join(data_path, "spatial")):
        loaded_kind = "visium-dir"
        info("reading Visium directory " + data_path)
        adata = sc.read_visium(data_path, var_names="gene_symbols")
        platform_info = {"platform": "visium", "input_kind": "visium-dir"}
    else:
        raise ValueError(
            "--input is a directory with no spatial/ folder; expected a Space Ranger "
            "output (spatial/ + matrix) for Visium, or pass a Stereo-seq .gef/.gem file"
        )
elif lower.endswith(".h5ad"):
    loaded_kind = "anndata"
    adata = sc.read_h5ad(data_path)
    if "spatial" not in adata.obsm and {"x", "y"} <= set(adata.obs):
        adata.obsm["spatial"] = adata.obs[["x", "y"]].to_numpy(dtype=float)
    platform_info = {"platform": "anndata", "input_kind": "anndata"}
elif lower.endswith(".gef"):
    loaded_kind = "stereoseq-gef"
    adata, platform_info = _read_gef_square(data_path, bin_size)
elif lower.endswith((".gem", ".gem.gz")):
    loaded_kind = "stereoseq-gem"
    adata, platform_info = _read_gem(data_path, bin_size, cellbin)
elif lower.endswith(".h5"):
    loaded_kind = "visium-h5"
    adata = sc.read_10x_h5(data_path, var_names="gene_symbols")
    adata.var_names_make_unique()
    platform_info = {"platform": "visium", "input_kind": "visium-h5"}
    info("read a .h5 matrix; if it has no tissue positions, pass the Space Ranger directory instead")
elif lower.endswith((".csv", ".tsv", ".csv.gz", ".tsv.gz")):
    loaded_kind = "spatial-table"
    sep = "," if lower.endswith((".csv", ".csv.gz")) else "\t"
    frame = pd.read_csv(data_path, sep=sep)
    bx, by = _bin_xy(frame["x"], frame["y"], bin_size)
    adata = _bins_to_adata(bx, by, frame["gene"].astype(str).to_numpy(),
                           None, frame["count"].astype(np.float64).to_numpy())
    platform_info = {"platform": "visium-hd", "input_kind": "spatial-table", "bin_size": bin_size}
else:
    raise ValueError(
        "unsupported --input; expected a Visium directory/.h5, a Stereo-seq .gef/.gem, "
        "a Visium HD bin table (.csv/.tsv with gene,x,y,count), or an .h5ad with spatial coordinates"
    )

adata.var_names_make_unique()
if "spatial" not in adata.obsm:
    raise ValueError(
        "loaded object has no obsm['spatial']; spatial statistics need coordinates. "
        "For Visium pass the Space Ranger directory (spatial/tissue_positions.csv)."
    )
if "counts" not in adata.layers:
    info("no counts layer; snapshotting X as counts")
    adata.layers["counts"] = adata.X.copy()
adata.uns["spatial_platform"] = platform_info
bin_info = {k: platform_info.get(k) for k in ("bin_size", "bin_auto", "available_bins", "cellbin")}

# HARD GATE: if min_genes would remove every spatial unit, say so here instead of
# failing later with an empty-object error.
gate = 0 if min_genes is None else min_genes
detected = np.asarray((adata.layers["counts"] > 0).sum(axis=1)).ravel()
print(f"genes detected per unit: min {detected.min()}, median {int(np.median(detected))}, max {detected.max()}")
if (detected > gate).sum() == 0:
    raise ValueError(
        f"No spatial unit has more than {gate} detected genes, so filtering would remove "
        "every unit. Lower min_genes in the config cell, check the bin size, or verify the "
        "matrix rows are genes."
    )
print(adata)