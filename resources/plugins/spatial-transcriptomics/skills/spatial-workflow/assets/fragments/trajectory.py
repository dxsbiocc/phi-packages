# Spatial trajectory (PAGA + diffusion pseudotime) and gene-program spatial
# scoring map developmental or state-transition gradients back onto tissue
# coordinates.
#
# 1. PAGA (Partition-based Graph Abstraction) quantifies expression-manifold
#    connectivity between clusters/cell types; diffusion pseudotime (DPT) orders
#    units along that manifold from a root state and projects pseudotime onto
#    obsm[coord_key].
# 2. Gene-program scoring (sc.tl.score_genes) evaluates user-supplied pathway or
#    signature gene sets on the unscaled lognorm layer and maps them in space.
import json
import matplotlib.pyplot as plt

if group_by not in adata.obs.columns:
    raise ValueError(f"group column '{group_by}' not in obs")
if not isinstance(adata.obs[group_by].dtype, pd.CategoricalDtype):
    adata.obs[group_by] = adata.obs[group_by].astype("category")

# --- Part 1: Gene-program spatial scoring (optional via gene_sets_path) ---
if gene_sets_path:
    if not os.path.exists(gene_sets_path):
        raise ValueError("gene_sets_path not found: " + gene_sets_path)
    with open(gene_sets_path, encoding="utf-8") as handle:
        gene_sets = json.load(handle)
    if not isinstance(gene_sets, dict) or not gene_sets:
        raise ValueError("gene_sets_path must be a non-empty JSON object mapping program -> [genes]")

    expr = adata.layers["lognorm"] if "lognorm" in adata.layers else adata.X
    score_ad = sc.AnnData(X=expr, obs=adata.obs[[group_by]].copy(), var=adata.var.copy())
    scored_programs = {}
    score_cols = []
    for prog_name, genes in gene_sets.items():
        present = [g for g in genes if g in score_ad.var_names]
        if len(present) < 2:
            analysis_notes.append(f"gene-set '{prog_name}' skipped: fewer than 2 genes present")
            continue
        col = f"score_{prog_name}"
        n_bins = min(25, max(2, score_ad.n_vars // 5))
        sc.tl.score_genes(
            score_ad,
            gene_list=present,
            score_name=col,
            ctrl_size=min(50, max(5, score_ad.n_vars // 4)),
            n_bins=n_bins,
            random_state=seed,
            use_raw=False,
        )
        adata.obs[col] = score_ad.obs[col].to_numpy(dtype=float)
        score_cols.append(col)
        by_grp = adata.obs.groupby(group_by, observed=True)[col].mean()
        scored_programs[str(prog_name)] = {
            "column": col,
            "genes_present": len(present),
            "genes_total": len(genes),
            "mean_by_group": {str(k): round(float(v), 4) for k, v in by_grp.items()},
        }

    score_results = {
        "gene_sets_path": os.path.basename(gene_sets_path),
        "programs_scored": len(scored_programs),
        "programs": scored_programs,
    }
    if score_cols and coord_key in adata.obsm:
        coords = np.asarray(adata.obsm[coord_key])
        n_show = min(4, len(score_cols))
        fig, axes = plt.subplots(1, n_show, figsize=(4.2 * n_show, 4.0), squeeze=False)
        for idx, col in enumerate(score_cols[:n_show]):
            ax = axes[0, idx]
            sc_plot = ax.scatter(coords[:, 0], coords[:, 1], c=adata.obs[col].to_numpy(), s=12, cmap="RdBu_r")
            ax.invert_yaxis()
            ax.set_aspect("equal", adjustable="datalim")
            ax.set_title(col)
            fig.colorbar(sc_plot, ax=ax, fraction=0.046, pad=0.04)
        fig.tight_layout()
        show_plot(fig, "trajectory-gene-scores")
else:
    score_results = {"skipped": "no gene_sets_path provided"}

# --- Part 2: PAGA + Diffusion Pseudotime mapped to tissue coordinates ---
categories = list(adata.obs[group_by].cat.categories)
if len(categories) < 2 or adata.n_obs < 10:
    info("trajectory needs >= 2 groups and >= 10 units; skipping PAGA/DPT")
    trajectory_results = {"skipped": "fewer than 2 groups or < 10 units"}
    analysis_notes.append("trajectory skipped: fewer than 2 groups or < 10 units")
else:
    # Ensure expression neighbors exist for PAGA and diffusion map.
    if "neighbors" not in adata.uns:
        if "X_pca" not in adata.obsm:
            n_pcs_traj = min(int(n_pcs), max(2, min(adata.n_obs, adata.n_vars) - 1))
            sc.pp.pca(adata, n_comps=n_pcs_traj, random_state=seed)
        sc.pp.neighbors(adata, n_neighbors=min(15, max(3, adata.n_obs - 1)), random_state=seed)

    import scipy.sparse as sp

    paga_disconnected = False
    try:
        sc.tl.paga(adata, groups=group_by)
    except ValueError as err:
        if "mismatching number of index arrays" not in str(err):
            raise
        # Known python-igraph + scipy>=1.16 bug when inter-cluster edge list is empty
        # (completely disconnected clusters in the expression kNN graph).
        paga_disconnected = True
        n_cats = len(categories)
        zero_mat = sp.csr_matrix((n_cats, n_cats), dtype=float)
        adata.uns["paga"] = {
            "groups": group_by,
            "connectivities": zero_mat,
            "connectivities_tree": zero_mat,
        }
        analysis_notes.append(
            "PAGA found no inter-cluster edges in the expression kNN graph (clusters are disconnected); "
            "DPT pseudotime may be infinite for disconnected components"
        )
    n_dcs = min(10, max(3, adata.n_obs - 2))
    sc.tl.diffmap(adata, n_comps=n_dcs)

    # Choose root cell: explicit root_cluster if given, otherwise first category with a note.
    if root_cluster is not None and str(root_cluster) in [str(c) for c in categories]:
        chosen_root = str(root_cluster)
        root_selection = "explicit"
    else:
        if root_cluster is not None:
            analysis_notes.append(f"root_cluster '{root_cluster}' not in {group_by}; fell back to auto root")
        chosen_root = str(categories[0])
        root_selection = "auto"
        analysis_notes.append(
            f"trajectory used auto root_cluster='{chosen_root}'; pass root_cluster to set biological direction"
        )

    root_mask = (adata.obs[group_by].astype(str) == chosen_root).to_numpy()
    root_indices = np.flatnonzero(root_mask)
    # Pick the cell in root_cluster closest to the cluster centroid in diffusion space.
    diffmap = np.asarray(adata.obsm["X_diffmap"])
    centroid = diffmap[root_indices].mean(axis=0)
    dists = np.linalg.norm(diffmap[root_indices] - centroid, axis=1)
    iroot = int(root_indices[int(np.argmin(dists))])
    adata.uns["iroot"] = iroot

    sc.tl.dpt(adata, n_dcs=min(n_dcs, diffmap.shape[1]))
    analysis_steps.append("trajectory")

    conn = adata.uns["paga"]["connectivities"]
    if hasattr(conn, "toarray"):
        conn = conn.toarray()
    conn = np.asarray(conn, dtype=float)
    paga_edges = []
    for i in range(len(categories)):
        for j in range(i + 1, len(categories)):
            weight = float(conn[i, j])
            if weight >= float(paga_threshold):
                paga_edges.append(
                    {"from": str(categories[i]), "to": str(categories[j]), "connectivity": round(weight, 4)}
                )
    paga_edges.sort(key=lambda e: -e["connectivity"])

    pt = adata.obs["dpt_pseudotime"].to_numpy(dtype=float)
    finite_pt = pt[np.isfinite(pt)]
    median_by_group = (
        adata.obs.loc[np.isfinite(pt)]
        .groupby(group_by, observed=True)["dpt_pseudotime"]
        .median()
    )
    trajectory_results = {
        "group_by": group_by,
        "root_cluster": chosen_root,
        "root_selection": root_selection,
        "iroot": iroot,
        "paga_threshold": float(paga_threshold),
        "n_paga_edges": len(paga_edges),
        "top_paga_edges": paga_edges[:10],
        "pseudotime_quantiles": {
            "p10": round(float(np.percentile(finite_pt, 10)), 4) if len(finite_pt) else None,
            "p50": round(float(np.percentile(finite_pt, 50)), 4) if len(finite_pt) else None,
            "p90": round(float(np.percentile(finite_pt, 90)), 4) if len(finite_pt) else None,
        },
        "median_pseudotime_by_group": {
            str(k): round(float(v), 4) for k, v in median_by_group.sort_values().items()
        },
    }
    print(f"trajectory (root={chosen_root} [{root_selection}]): {len(paga_edges)} PAGA edges >= {paga_threshold}")
    for grp, val in trajectory_results["median_pseudotime_by_group"].items():
        print(f"  {grp}: median DPT {val}")

    # Plot PAGA graph.
    try:
        sc.pl.paga(adata, threshold=float(paga_threshold), show=False)
        show_plot(plt.gcf(), "trajectory-paga")
    except Exception as error:
        info(f"PAGA plot skipped: {error}")

    # Plot diffusion pseudotime mapped onto spatial coordinates.
    if coord_key in adata.obsm:
        coords = np.asarray(adata.obsm[coord_key])
        fig, ax = plt.subplots(figsize=(5.5, 5.0))
        valid = np.isfinite(pt)
        sc_plot = ax.scatter(coords[valid, 0], coords[valid, 1], c=pt[valid], s=14, cmap="viridis")
        ax.invert_yaxis()
        ax.set_aspect("equal", adjustable="datalim")
        ax.set_title(f"Spatial diffusion pseudotime (root: {chosen_root})")
        fig.colorbar(sc_plot, ax=ax, fraction=0.046, pad=0.04, label="dpt_pseudotime")
        show_plot(fig, "trajectory-spatial-dpt")
