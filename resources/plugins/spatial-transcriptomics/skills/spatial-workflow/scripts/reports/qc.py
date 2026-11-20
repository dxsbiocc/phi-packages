# JSON contract for st_py_qc. Reads the bookkeeping the fragments left behind.
before = {"units": obj_before.n_obs, "genes": obj_before.n_vars, "per_sample": None}
if sample_column in obj_before.obs.columns:
    before["per_sample"] = [
        {"sample": str(name), "units": int(count)}
        for name, count in obj_before.obs[sample_column].value_counts().items()
    ]

metric_cols = [c for c in ("n_genes_by_counts", "total_counts", "pct_counts_mt") if c in obj_before.obs]
distribution = []
for column in metric_cols:
    values = obj_before.obs[column]
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
    "units": adata.n_obs,
    "genes": adata.n_vars,
    "removed_units": before["units"] - adata.n_obs,
    "removed_genes": before["genes"] - adata.n_vars,
    "kept_fraction": round(adata.n_obs / max(1, before["units"]), 4),
}

warnings = []
if adata.n_obs == 0:
    warnings.append("all spatial units were removed; relax the thresholds or check the bin size")
if 0 < after["kept_fraction"] < 0.5:
    warnings.append(
        f"only {after['kept_fraction'] * 100:.1f}% of units survived; "
        "check whether a threshold is too strict or the bin is too small"
    )
if "pct_counts_mt" not in metric_cols:
    warnings.append(
        f"no gene matched mt_pattern ({mt_pattern}); the mitochondrial gate was skipped "
        "(use ^mt- for mouse)"
    )

emit(
    {
        "input": tool_input,
        "output": tool_output,
        "platform": platform_info,
        "bin": bin_info,
        "before": before,
        "after": after,
        "distribution": distribution,
        "filters": qc_filters,
        "features_scored": qc_scored,
        "figures": [
            f
            for f in (
                figure_paths.get("qc-violin-before"),
                figure_paths.get("qc-violin-after"),
                figure_paths.get("qc-spatial-before"),
                figure_paths.get("qc-spatial-after"),
            )
            if f
        ],
        "warnings": warnings,
    }
)