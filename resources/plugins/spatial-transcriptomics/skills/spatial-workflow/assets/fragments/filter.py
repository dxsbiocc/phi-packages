# Thresholds come from the config cell. Read the describe() output and the
# violin panels above, then edit those values rather than accepting them blindly:
#   n_genes_by_counts > min_genes    removes empty / off-tissue units
#   total_counts      > min_counts   removes empty units
#   total_counts      < max_counts   removes oversized units (e.g. merged bins)
#   pct_counts_mt     < mt_max       removes damaged tissue
# A None bound skips that gate. Each gate is applied in sequence and reports how
# many units IT removed, so a surprising total can be traced to one threshold.
# Spatial units are not droplets: there is no doublet gate here.
obj_before = adata
qc_filters = []
kept = np.ones(adata.n_obs, dtype=bool)


def apply_filter(name, condition):
    global kept
    condition = np.asarray(condition, dtype=bool)
    qc_filters.append({"filter": name, "removed": int((~condition).sum())})
    kept &= condition
    print(f"{name:<28} removed {int((~condition).sum())} units")


if min_genes is not None:
    apply_filter(f"n_genes_by_counts >= {min_genes}", adata.obs["n_genes_by_counts"] >= min_genes)
if min_counts is not None:
    apply_filter(f"total_counts >= {min_counts}", adata.obs["total_counts"] >= min_counts)
if max_counts is not None:
    apply_filter(f"total_counts <= {max_counts}", adata.obs["total_counts"] <= max_counts)
if mt_max is not None:
    if "pct_counts_mt" in adata.obs:
        apply_filter(f"pct_counts_mt <= {mt_max}", adata.obs["pct_counts_mt"] <= mt_max)
    else:
        print(f"mt_max={mt_max} skipped: no gene matched mt_pattern ({mt_pattern}); "
              "use ^mt- for mouse or check the identifier type")

adata = adata[kept, :].copy()
print(f"kept {adata.n_obs} of {obj_before.n_obs} units "
      f"({100 * adata.n_obs / max(1, obj_before.n_obs):.1f}%)")

# Keeping under half the units usually means a threshold is too strict for this
# dataset, or the bin size is too small. Re-read the distributions before continuing.
if adata.n_obs == 0:
    print("WARNING: every unit was filtered out; relax the thresholds.")
elif adata.n_obs < obj_before.n_obs / 2:
    print("WARNING: more than half the units were removed; reconsider the thresholds in the config cell.")