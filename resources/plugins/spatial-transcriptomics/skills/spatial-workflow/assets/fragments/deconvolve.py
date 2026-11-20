# Deconvolution estimates the cell-type proportions of each spatial unit. A Visium
# spot or a Stereo-seq bin mixes several cells, so its expression is a weighted sum
# of cell-type profiles. With a single-cell reference we can solve for those weights.
#
# This is a non-negative least-squares (NNLS) estimate, NOT a trained model: for each
# unit it finds non-negative weights w minimising ||S w - y||, where S is the
# reference signature (mean normalized profile per cell type) and y the unit's
# normalized expression, then rescales w to sum to 1. It is fast, dependency-free
# (scipy only) and honest about being a linear approximation. For model-based
# deconvolution (cell2location, tangram-sc, stereopy) use a project environment —
# those are pip-only and cannot enter this pure-conda locked environment.
from scipy import sparse
from scipy.optimize import nnls

if not reference_path:
    info("no reference_path; skipping deconvolution. Pass a single-cell .h5ad to estimate proportions.")
    deconv_results = {"skipped": "no reference_path provided"}
    analysis_notes.append("deconvolution skipped: no reference")
else:
    if not os.path.exists(reference_path):
        raise ValueError("reference not found: " + reference_path)
    ref = sc.read_h5ad(reference_path)
    if reference_cell_type_key not in ref.obs.columns:
        raise ValueError(
            f"reference has no obs column '{reference_cell_type_key}'; available: "
            + ", ".join(ref.obs.columns)
        )
    if deconv_method != "nnls":
        raise ValueError("deconv_method must be nnls in-plugin; cell2location/tangram need a project environment")

    # Work in linear, library-size-normalized space. Prefer the reference counts
    # layer; without it, assume X is log1p-normalized (the scanpy default) and invert.
    ref_expr = ref.layers["counts"] if "counts" in ref.layers else ref.X
    if not sparse.issparse(ref_expr):
        ref_expr = sparse.csr_matrix(ref_expr)
    ref_expr = ref_expr.tocsr()
    if "counts" not in ref.layers:
        ref_expr.data = np.expm1(ref_expr.data)
        analysis_notes.append("reference had no counts layer; assumed log1p X and inverted with expm1")
    ref_lib = np.asarray(ref_expr.sum(axis=1)).ravel()
    ref_lib[ref_lib == 0] = 1.0
    ref_norm = sparse.diags(1.0 / ref_lib) @ ref_expr

    # Signature genes: shared with the reference, preferring spatial HVGs, ranked by
    # spatial variance and capped so the per-unit NNLS stays tractable.
    shared = [g for g in adata.var_names if g in set(ref.var_names)]
    if not shared:
        raise ValueError("no shared genes between spatial object and reference; check identifier type")
    if "highly_variable" in adata.var:
        hvg = set(adata.var_names[adata.var["highly_variable"]])
        candidates = [g for g in shared if g in hvg] or shared
    else:
        candidates = shared
    counts = adata.layers["counts"]
    if not sparse.issparse(counts):
        counts = sparse.csr_matrix(counts)
    g_mean = np.asarray(counts.mean(axis=0)).ravel()
    g_var = np.asarray(counts.multiply(counts).mean(axis=0)).ravel() - g_mean**2
    sp_idx = {g: i for i, g in enumerate(adata.var_names)}
    candidates = sorted(candidates, key=lambda g: -g_var[sp_idx[g]])
    sig_genes = candidates[:500]

    ref_idx = {g: i for i, g in enumerate(ref.var_names)}
    sig_ref = [ref_idx[g] for g in sig_genes]
    types = sorted_levels(ref.obs[reference_cell_type_key].astype(str))
    signature = np.zeros((len(sig_genes), len(types)))
    labels = ref.obs[reference_cell_type_key].astype(str).to_numpy()
    for j, cell_type in enumerate(types):
        mask = labels == cell_type
        signature[:, j] = np.asarray(ref_norm[mask][:, sig_ref].mean(axis=0)).ravel()

    lib = np.asarray(counts.sum(axis=1)).ravel()
    lib[lib == 0] = 1.0
    norm = sparse.diags(1.0 / lib) @ counts
    sig_sp = [sp_idx[g] for g in sig_genes]
    Y = np.asarray(norm[:, sig_sp].todense())

    proportions = np.zeros((adata.n_obs, len(types)))
    for i in range(adata.n_obs):
        w, _ = nnls(signature, Y[i])
        total = w.sum()
        proportions[i] = w / total if total > 0 else 0.0

    adata.obsm["deconvolution"] = proportions
    for j, cell_type in enumerate(types):
        adata.obs[f"{cell_type}_prop"] = proportions[:, j]
    analysis_steps.append("deconvolve")

    mean_prop = proportions.mean(axis=0)
    deconv_results = {
        "method": "nnls",
        "reference": os.path.basename(reference_path),
        "cell_types": types,
        "signature_genes": len(sig_genes),
        "n_units": int(adata.n_obs),
        "mean_proportions": {ct: round(float(mean_prop[j]), 4) for j, ct in enumerate(types)},
        "note": "linear NNLS estimate, not a validated model; use cell2location in a project env for model-based deconvolution",
    }
    print(f"deconvolved {adata.n_obs} units into {len(types)} cell types over {len(sig_genes)} signature genes")
    print("mean proportions:", deconv_results["mean_proportions"])

    top_types = [types[j] for j in np.argsort(-mean_prop)[:3]]
    if coord_key in adata.obsm:
        show_plot(
            sc.pl.embedding(adata, basis=coord_key, color=[f"{ct}_prop" for ct in top_types], show=False),
            "deconv-proportions-spatial",
        )