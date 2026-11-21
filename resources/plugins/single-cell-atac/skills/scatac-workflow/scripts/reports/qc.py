# JSON contract for scatac_qc. Reads the bookkeeping the fragments left behind.
# obj_before is a plain-number snapshot (snapatac2 filters the backed object in place).
before = obj_before or {"cells": 0, "genes": 0, "per_sample": None}


def _dist(values):
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    if not len(v):
        return None
    return {
        "min": round(float(v.min()), 4),
        "q01": round(float(np.quantile(v, quantiles[0])), 4),
        "median": round(float(np.median(v)), 4),
        "q99": round(float(np.quantile(v, quantiles[-1])), 4),
        "max": round(float(v.max()), 4),
    }


distribution = []
for column in ("n_fragment", "frac_dup", "frac_mito", "tsse"):
    if obs_has(adata, column):
        d = _dist(obs_col(adata, column))
        if d:
            distribution.append({"metric": column, **d})

after = {
    "cells": adata.n_obs,
    "genes": adata.n_vars,
    "removed_cells": before["cells"] - adata.n_obs,
    "kept_fraction": round(adata.n_obs / max(1, before["cells"]), 4),
}

warnings = []
if adata.n_obs == 0:
    warnings.append("all cells were removed; relax min_counts/min_tsse or check the genome build")
elif 0 < after["kept_fraction"] < 0.5:
    warnings.append(
        f"only {after['kept_fraction'] * 100:.1f}% of cells survived; check whether a threshold is too strict"
    )
if "tsse" not in [d["metric"] for d in distribution]:
    warnings.append("no tsse computed; pass gene_anno to enable TSS enrichment, the best ATAC signal-to-noise metric")

emit(
    {
        "input": tool_input,
        "output": tool_output,
        "genome": genome_info,
        "before": before,
        "after": after,
        "distribution": distribution,
        "filters": qc_filters,
        "features_scored": qc_scored,
        "figures": [
            f
            for f in (
                figure_paths.get("qc-frag-size-before"),
                figure_paths.get("qc-frag-size-after"),
                figure_paths.get("qc-tsse-before"),
                figure_paths.get("qc-tsse-after"),
                figure_paths.get("qc-metrics-before"),
                figure_paths.get("qc-metrics-after"),
            )
            if f
        ],
        "warnings": warnings,
        "notes": analysis_notes,
    }
)