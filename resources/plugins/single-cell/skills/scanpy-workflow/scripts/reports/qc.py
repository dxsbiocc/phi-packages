# JSON contract for sc_py_qc. Reads the bookkeeping the fragments left behind.
before = {"cells": obj_before.n_obs, "genes": obj_before.n_vars, "per_sample": None}
if sample_column in obj_before.obs.columns:
    before["per_sample"] = [
        {"sample": str(name), "cells": int(count)}
        for name, count in obj_before.obs[sample_column].value_counts().items()
    ]

metrics = obj_before.obs[["n_genes_by_counts", "total_counts", "pct_counts_mt"]]
distribution = []
for column in metrics.columns:
    values = metrics[column]
    distribution.append(
        {
            "metric": column,
            "min": round(float(values.min()), 3),
            "q01": round(float(values.quantile(quantiles[0])), 3),
            "median": round(float(values.median()), 3),
            "q99": round(float(values.quantile(quantiles[-1])), 3),
            "max": round(float(values.max()), 3),
        }
    )

after = {
    "cells": adata.n_obs,
    "genes": adata.n_vars,
    "removed_cells": before["cells"] - adata.n_obs,
    "removed_genes": before["genes"] - adata.n_vars,
    "kept_fraction": round(adata.n_obs / max(1, before["cells"]), 4),
}

warnings = []
if adata.n_obs == 0:
    warnings.append("all cells were removed; relax the thresholds")
if 0 < after["kept_fraction"] < 0.5:
    warnings.append(
        f"only {after['kept_fraction'] * 100:.1f}% of cells survived; "
        "check whether a threshold is too strict"
    )
if detect_doublets and doublet_result and doublet_result.get("error"):
    warnings.append("doublet detection unavailable: " + doublet_result["error"])

emit(
    {
        "input": tool_input,
        "output": tool_output,
        "before": before,
        "after": after,
        "distribution": distribution,
        "filters": qc_filters,
        "features_scored": qc_scored,
        "doublets": doublet_result,
        "figures": [
            f
            for f in (figure_paths.get("qc-violin-before"), figure_paths.get("qc-violin-after"))
            if f
        ],
        "warnings": warnings,
    }
)