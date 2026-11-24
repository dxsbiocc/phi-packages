# Spatial niche (tissue domain) segmentation partitions tissue into recurring
# microenvironments by aggregating each unit's spatial neighbourhood — a
# lightweight, pure-conda formulation of UTAG (message passing over the spatial
# graph followed by Leiden clustering).
#
# Two modes are supported:
# - composition (default): aggregate the local cell-type/cluster composition
#   around each unit (or deconvolution proportions when niche_use_deconvolution
#   is true and obsm['deconvolution'] exists), then cluster units with similar
#   neighbourhood composition into niches.
# - expression: smooth the lognorm HVG expression over the spatial graph
#   (W @ X) and cluster the smoothed representation into spatial domains.
#
# When group_column names a condition/stage column with >= 2 conditions, we also
# compute per-condition niche proportions and cross-stage shifts.
import scipy.sparse as sp
import matplotlib.pyplot as plt
import seaborn as sns

if "spatial_connectivities" not in adata.obsp:
    raise ValueError("no spatial graph; run the spatial-neighbors step first")

_niche_col = niche_column or "niche"
_niche_mode = (niche_mode or "composition").lower()
if _niche_mode not in ("composition", "expression"):
    raise ValueError(f"unsupported niche_mode '{_niche_mode}'; use 'composition' or 'expression'")

# Build row-normalized (A + I)^hops spatial aggregation operator.
A = adata.obsp["spatial_connectivities"]
if not sp.issparse(A):
    A = sp.csr_matrix(A)
A = A.astype(float)
W = A + sp.eye(adata.n_obs, format="csr", dtype=float)
row_sums = np.asarray(W.sum(axis=1)).ravel()
row_sums[row_sums == 0] = 1.0
W = sp.diags(1.0 / row_sums) @ W
_hops = max(1, int(niche_hops))
if _hops > 1:
    W_hop = W.copy()
    for _ in range(_hops - 1):
        W_hop = W_hop @ W
    W = W_hop

feature_names = None
if _niche_mode == "composition":
    if niche_use_deconvolution and "deconvolution" in adata.obsm:
        deconv = adata.obsm["deconvolution"]
        if isinstance(deconv, pd.DataFrame):
            feat_mat = deconv.to_numpy(dtype=float)
            feature_names = [str(c) for c in deconv.columns]
        else:
            feat_mat = np.asarray(deconv, dtype=float)
            feature_names = [f"type_{i}" for i in range(feat_mat.shape[1])]
        comp_source = "deconvolution"
    else:
        if group_by not in adata.obs.columns:
            raise ValueError(f"group column '{group_by}' not in obs")
        dummies = pd.get_dummies(adata.obs[group_by].astype(str), dtype=float)
        feat_mat = dummies.to_numpy(dtype=float)
        feature_names = [str(c) for c in dummies.columns]
        comp_source = f"obs:{group_by}"
    agg_mat = np.asarray(W @ feat_mat, dtype=float)
    adata.obsm["X_niche"] = agg_mat
    niche_ad = sc.AnnData(X=agg_mat, obs=adata.obs.copy())
    n_pcs_niche = min(15, max(2, agg_mat.shape[1] - 1), max(2, adata.n_obs - 1))
    if agg_mat.shape[1] > 3 and adata.n_obs > n_pcs_niche + 1:
        sc.pp.pca(niche_ad, n_comps=n_pcs_niche, random_state=seed)
        sc.pp.neighbors(niche_ad, n_neighbors=min(15, max(3, adata.n_obs - 1)), use_rep="X_pca", random_state=seed)
    else:
        sc.pp.neighbors(niche_ad, n_neighbors=min(15, max(3, adata.n_obs - 1)), use_rep="X", random_state=seed)
else:
    expr = adata.layers["lognorm"] if "lognorm" in adata.layers else adata.X
    if "highly_variable" in adata.var.columns and int(adata.var["highly_variable"].sum()) >= 10:
        hvg_mask = adata.var["highly_variable"].to_numpy()
        expr_sub = expr[:, hvg_mask]
    else:
        expr_sub = expr
    agg_expr = W @ expr_sub
    if sp.issparse(agg_expr):
        agg_expr = agg_expr.toarray()
    agg_expr = np.asarray(agg_expr, dtype=float)
    niche_ad = sc.AnnData(X=agg_expr, obs=adata.obs.copy())
    n_pcs_niche = min(int(n_pcs), max(2, agg_expr.shape[1] - 1), max(2, adata.n_obs - 1))
    sc.pp.scale(niche_ad, max_value=10)
    sc.pp.pca(niche_ad, n_comps=n_pcs_niche, random_state=seed)
    adata.obsm["X_niche"] = niche_ad.obsm["X_pca"].copy()
    sc.pp.neighbors(niche_ad, n_neighbors=min(15, max(3, adata.n_obs - 1)), use_rep="X_pca", random_state=seed)
    comp_source = "lognorm_smoothed"

sc.tl.leiden(
    niche_ad,
    resolution=float(niche_resolution),
    key_added=_niche_col,
    random_state=seed,
    flavor="igraph",
    n_iterations=2,
    directed=False,
)
adata.obs[_niche_col] = pd.Categorical([f"niche_{c}" for c in niche_ad.obs[_niche_col].astype(str)])
analysis_steps.append("niche")

# Summarize niche sizes and constituent group_by composition.
niche_counts = adata.obs[_niche_col].value_counts().sort_index()
niche_composition = {}
if group_by in adata.obs.columns:
    ct = pd.crosstab(adata.obs[_niche_col], adata.obs[group_by].astype(str), normalize="index")
    for n_name, row in ct.iterrows():
        top_members = row.sort_values(ascending=False).head(5)
        niche_composition[str(n_name)] = {str(k): round(float(v), 4) for k, v in top_members.items() if v > 0}
else:
    ct = None

# Cross-condition / stage niche dynamics when group_column is provided.
condition_dynamics = None
if group_column and group_column in adata.obs.columns:
    conds = sorted_levels(adata.obs[group_column].astype(str))
    if len(conds) >= 2:
        dyn_ct = pd.crosstab(adata.obs[group_column].astype(str), adata.obs[_niche_col], normalize="index")
        dyn_ct = dyn_ct.loc[[c for c in conds if c in dyn_ct.index]]
        shifts = {}
        for col in dyn_ct.columns:
            vals = dyn_ct[col]
            shifts[str(col)] = {
                "by_condition": {str(c): round(float(vals.loc[c]), 4) for c in dyn_ct.index},
                "range": round(float(vals.max() - vals.min()), 4),
            }
        condition_dynamics = {
            "group_column": group_column,
            "conditions": [str(c) for c in dyn_ct.index],
            "niche_proportions": shifts,
        }
        fig, ax = plt.subplots(figsize=(max(4.5, len(dyn_ct.index) * 1.4), 4.0))
        dyn_ct.plot(kind="bar", stacked=True, ax=ax, colormap="tab20")
        ax.set_ylabel("Niche proportion")
        ax.set_xlabel(group_column)
        ax.set_title(f"Tissue niche dynamics across {group_column}")
        ax.legend(title=_niche_col, bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
        show_plot(fig, "niche-condition-dynamics")

niche_results = {
    "niche_column": _niche_col,
    "mode": _niche_mode,
    "source": comp_source,
    "hops": _hops,
    "resolution": float(niche_resolution),
    "n_niches": int(len(niche_counts)),
    "niche_sizes": {str(k): int(v) for k, v in niche_counts.items()},
    "composition": niche_composition,
    "condition_dynamics": condition_dynamics,
}
adata.uns["niche_summary"] = {
    "niche_column": _niche_col,
    "mode": _niche_mode,
    "hops": _hops,
    "resolution": float(niche_resolution),
}

print(f"niche ({_niche_mode}, hops={_hops}, res={niche_resolution}): {len(niche_counts)} niches")
for k, v in niche_counts.items():
    print(f"  {k}: {int(v)} units")

# Spatial niche map.
coords = np.asarray(adata.obsm[coord_key])
fig, ax = plt.subplots(figsize=(5.5, 5.0))
niche_cats = list(adata.obs[_niche_col].cat.categories)
palette = plt.cm.tab20(np.linspace(0, 1, max(1, len(niche_cats))))
for idx, cat in enumerate(niche_cats):
    mask = (adata.obs[_niche_col] == cat).to_numpy()
    ax.scatter(coords[mask, 0], coords[mask, 1], s=14, color=palette[idx], label=str(cat), alpha=0.9)
ax.invert_yaxis()
ax.set_aspect("equal", adjustable="datalim")
ax.set_title(f"Spatial niches ({_niche_col})")
ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False, markerscale=1.4)
show_plot(fig, "niche-spatial")

# Niche composition heatmap.
if ct is not None and ct.shape[0] >= 1 and ct.shape[1] >= 2:
    fig, ax = plt.subplots(figsize=(max(4.5, ct.shape[1] * 0.7), max(3.5, ct.shape[0] * 0.5)))
    sns.heatmap(ct, annot=True, fmt=".2f", cmap="YlGnBu", ax=ax, cbar_kws={"label": "Proportion"})
    ax.set_title(f"Niche composition by {group_by}")
    show_plot(fig, "niche-composition")
