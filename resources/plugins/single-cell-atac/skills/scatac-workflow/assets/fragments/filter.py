# Threshold gates on the per-cell QC metrics. Read the qc-metrics distributions
# and the qc-plots panels first, then set these in the config cell — a threshold
# is a decision, not a default. A None bound disables that gate. Each gate reports
# how many cells fail IT, so a surprising total traces to one threshold.
# Doublet removal is NOT here: scrublet needs the selected feature matrix, so it
# runs in the features step. This keeps the QC module lightweight (no tile matrix).
qc_filters = []
_nfrag = obs_col(adata, "n_fragment")
_tsse = obs_col(adata, "tsse")
_before = adata.n_obs


def _gate(name, mask):
    if mask is None:
        return
    removed = int((~np.asarray(mask, dtype=bool)).sum())
    qc_filters.append({"filter": name, "removed": removed})
    print(f"{name:<28} fails {removed} cells")


if min_counts is not None and _nfrag is not None:
    _gate(f"n_fragment >= {min_counts}", _nfrag >= min_counts)
if max_counts is not None and _nfrag is not None:
    _gate(f"n_fragment <= {max_counts}", _nfrag <= max_counts)

_eff_min_tsse = min_tsse
if min_tsse is not None:
    if _tsse is None:
        analysis_notes.append(f"min_tsse={min_tsse} skipped: tsse not computed (needs gene_anno)")
        info(f"min_tsse={min_tsse} skipped: no tsse column")
        _eff_min_tsse = None
    else:
        _gate(f"tsse >= {min_tsse}", _tsse >= min_tsse)
if max_tsse is not None and _tsse is not None:
    _gate(f"tsse <= {max_tsse}", _tsse <= max_tsse)

# None disables a gate (snapatac2 treats a falsy bound as off). filter_cells
# subsets the backed object in place; obj_before was snapshotted in qc-metrics.
snap.pp.filter_cells(
    adata, min_counts=min_counts, min_tsse=_eff_min_tsse,
    max_counts=max_counts, max_tsse=max_tsse, inplace=True,
)
print(f"kept {adata.n_obs} of {_before} cells ({100 * adata.n_obs / max(1, _before):.1f}%)")
if adata.n_obs == 0:
    raise ValueError(
        "every cell was filtered out; relax min_counts/min_tsse in the config cell, "
        "or check that the genome build matches the fragments."
    )
if adata.n_obs < _before / 2:
    analysis_notes.append(
        f"only {100 * adata.n_obs / max(1, _before):.1f}% of cells survived; reconsider the thresholds"
    )
    print("WARNING: more than half the cells were removed; reconsider the thresholds.")