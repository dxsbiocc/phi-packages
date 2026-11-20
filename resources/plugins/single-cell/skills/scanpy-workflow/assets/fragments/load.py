# data_path and project_name are set in the config cell above.
lower = data_path.lower()
if os.path.isdir(data_path):
    marker_files = ("features.tsv.gz", "genes.tsv.gz", "features.tsv", "genes.tsv")
    has_markers = any(os.path.exists(os.path.join(data_path, m)) for m in marker_files)
    has_matrix = os.path.exists(os.path.join(data_path, "matrix.mtx.gz")) or os.path.exists(
        os.path.join(data_path, "matrix.mtx")
    )
    nested = os.path.join(data_path, "filtered_feature_bc_matrix")
    if has_markers and has_matrix:
        loaded_kind = "10x-matrix-dir"
        directory = data_path
    elif os.path.isdir(nested):
        loaded_kind = "10x-matrix-dir"
        directory = nested
    else:
        raise ValueError(
            "--input is a directory with no 10x matrix files (matrix.mtx + features.tsv); "
            "expected a Cell Ranger filtered_feature_bc_matrix or raw_feature_bc_matrix folder"
        )
    info("reading 10x matrix directory " + directory)
    adata = sc.io.read_10x_mtx(directory, var_names="gene_symbols")
elif lower.endswith(".h5ad"):
    loaded_kind = "anndata"
    adata = sc.io.read_h5ad(data_path)
elif lower.endswith(".h5"):
    loaded_kind = "10x-h5"
    adata = sc.io.read_10x_h5(data_path)
elif lower.endswith(".loom"):
    loaded_kind = "loom"
    adata = sc.read_loom(data_path)
elif lower.endswith((".csv", ".tsv", ".txt")):
    loaded_kind = "dense-table"
    adata = sc.io.read_csv(data_path, delimiter="," if lower.endswith(".csv") else "\t")
elif lower.endswith((".mtx", ".mtx.gz")):
    raise ValueError(
        "--input " + data_path + " is a bare Matrix Market file without feature/barcode "
        "names; pass its 10x directory instead so cell and gene identifiers survive"
    )
else:
    raise ValueError(
        "unsupported --input; expected a .h5ad, a 10x matrix directory or .h5, "
        "a .loom, or a .csv/.tsv count table"
    )
adata.var_names_make_unique()

# Keep raw counts in a layer: QC, HVG (seurat_v3) and pseudobulk all need them.
if "counts" not in adata.layers:
    info("no counts layer; snapshotting X as counts")
    adata.layers["counts"] = adata.X.copy()

# HARD GATE: if min_genes would remove every cell, say so here instead of
# failing later with an empty-object error.
gate = 0 if min_genes is None else min_genes
detected = np.asarray((adata.layers["counts"] > 0).sum(axis=1)).ravel()
print(f"genes detected per cell: min {detected.min()}, median {int(np.median(detected))}, max {detected.max()}")
if (detected > gate).sum() == 0:
    raise ValueError(
        f"No cell has more than {gate} detected genes, so filtering would remove "
        "every cell. This usually means a filtered/subset matrix, a different feature "
        "ID type, or an unusually shallow library. Lower min_genes in the config cell "
        "(or --feature-low for the script tool) or check that the matrix rows are genes."
    )
print(adata)