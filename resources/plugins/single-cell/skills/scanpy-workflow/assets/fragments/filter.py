# Thresholds come from the config cell. Read the describe() output and the
# violin panels above, then edit those values rather than accepting them blindly:
#   n_genes_by_counts > min_genes    removes empty droplets
#   n_genes_by_counts < max_genes    removes likely doublets / oversized cells
#   total_counts      > min_counts   removes empty droplets
#   total_counts      < max_counts   removes likely doublets
#   pct_counts_mt     < mt_max       removes stressed cells
# A None bound skips that gate. Each gate is applied in sequence and reports how
# many cells IT removed, so a surprising total can be traced to one threshold.
obj_before = adata
qc_filters = []
kept = np.ones(adata.n_obs, dtype=bool)


def apply_filter(name, condition):
    global kept
    condition = np.asarray(condition, dtype=bool)
    qc_filters.append({"filter": name, "removed": int((~condition).sum())})
    kept &= condition
    print(f"{name:<28} removed {int((~condition).sum())} cells")


if min_genes is not None:
    apply_filter(f"n_genes_by_counts >= {min_genes}", adata.obs["n_genes_by_counts"] >= min_genes)
if max_genes is not None:
    apply_filter(f"n_genes_by_counts <= {max_genes}", adata.obs["n_genes_by_counts"] <= max_genes)
if min_counts is not None:
    apply_filter(f"total_counts >= {min_counts}", adata.obs["total_counts"] >= min_counts)
if max_counts is not None:
    apply_filter(f"total_counts <= {max_counts}", adata.obs["total_counts"] <= max_counts)
if mt_max is not None:
    apply_filter(f"pct_counts_mt <= {mt_max}", adata.obs["pct_counts_mt"] <= mt_max)
if detect_doublets and "predicted_doublet" in adata.obs.columns:
    apply_filter("predicted_doublet == False", ~adata.obs["predicted_doublet"].astype(bool))

adata = adata[kept, :].copy()
print(f"kept {adata.n_obs} of {obj_before.n_obs} cells "
      f"({100 * adata.n_obs / max(1, obj_before.n_obs):.1f}%)")

# Keeping under half the cells usually means a threshold is too strict for this
# dataset. Re-read the distributions before continuing.
if adata.n_obs == 0:
    print("WARNING: every cell was filtered out; relax the thresholds.")
elif adata.n_obs < obj_before.n_obs / 2:
    print("WARNING: more than half the cells were removed; reconsider the thresholds in the config cell.")