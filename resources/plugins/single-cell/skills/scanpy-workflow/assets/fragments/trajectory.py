# Pseudotime orders cells along a inferred progression; it is a relative rank,
# not a clock, and it is only as good as the root you choose. Pick the root from
# prior biology (stem/progenitor markers), not from the plot's left edge. PAGA
# gives the coarse cluster graph; diffusion pseudotime gives the per-cell order.
if group_by not in adata.obs.columns:
    raise ValueError(
        f"group column '{group_by}' not in obs; available: " + ", ".join(adata.obs.columns)
    )
if "neighbors" not in adata.uns:
    raise ValueError("no neighbour graph; run sc_py_preprocess before sc_py_trajectory")

levels = sorted_levels(adata.obs[group_by])
if root_cluster is None:
    root_cluster = levels[0]
elif root_cluster not in levels:
    raise ValueError(
        f"root cluster '{root_cluster}' not in {group_by}; available: " + ", ".join(levels)
    )

sc.tl.paga(adata, groups=group_by)
groups = adata.obs[group_by].astype(str)
root_cells = np.flatnonzero((groups == root_cluster).to_numpy())
if len(root_cells) == 0:
    raise ValueError(f"root cluster '{root_cluster}' has no cells")
adata.uns["iroot"] = int(root_cells[0])

sc.tl.diffmap(adata, n_comps=max(n_dcs, 15))
sc.tl.dpt(adata, n_dcs=n_dcs)

pseudotime = adata.obs["dpt_pseudotime"]
finite = pseudotime[np.isfinite(pseudotime)]
traj_info = {
    "group_by": group_by,
    "root_cluster": root_cluster,
    "n_dcs": n_dcs,
    "pseudotime": {
        "min": round(float(finite.min()), 4) if len(finite) else None,
        "max": round(float(finite.max()), 4) if len(finite) else None,
        "unreachable_cells": int(adata.n_obs - len(finite)),
    },
}
show_plot(sc.pl.paga(adata, show=False), "trajectory-paga")
if "X_umap" in adata.obsm:
    show_plot(
        sc.pl.umap(adata, color=[group_by, "dpt_pseudotime"], show=False),
        "trajectory-pseudotime",
    )
print(f"pseudotime rooted in cluster {root_cluster}; edit root_cluster if biology says another cluster starts the process")