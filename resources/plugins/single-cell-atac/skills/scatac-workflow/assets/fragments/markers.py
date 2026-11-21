# Marker regions are the peaks most specifically accessible in each cluster
# (marker_regions tests every cluster against the rest). Differential accessibility
# (diff_test) compares two clusters peak-by-peak with a logistic regression, giving
# log fold change and significance — the scATAC analogue of differential expression.
# Both run on the peak matrix (peak_adata) when peaks were called, otherwise on the
# tile matrix. diff_test takes cell INDICES, not cluster labels, so each cluster's
# cells are selected from the group column first.
from itertools import combinations

target = peak_adata if peak_adata is not None else adata
group = group_by if obs_has(target, group_by) else (cluster_columns[0] if cluster_columns else None)
if group is None:
    analysis_notes.append("markers skipped: no cluster column; run the cluster step first")
    info("markers skipped: no cluster column")
else:
    try:
        _mk = snap.tl.marker_regions(target, groupby=group, pvalue=marker_pvalue)
        markers = {str(k): list(v) for k, v in _mk.items()}
        os.makedirs(tables_dir, exist_ok=True)
        marker_tables = []
        for _cl, _regions in markers.items():
            _path = os.path.join(tables_dir, f"markers_{_cl}.txt")
            with open(_path, "w", encoding="utf-8") as handle:
                handle.write("\n".join(_regions) + ("\n" if _regions else ""))
            marker_tables.append({"cluster": _cl, "n_markers": len(_regions), "path": _path})
        analysis_steps.append("markers")
        print(f"marker regions per cluster: { {k: len(v) for k, v in markers.items()} }")

        # Pairwise diff_test is O(clusters^2) logistic regressions and is slow on a
        # real peak matrix, so it is opt-in (run_diff). marker_regions above already
        # gives each cluster's specific peaks against the rest; enable run_diff only
        # when a specific cluster-vs-cluster contrast is needed.
        _labels = obs_col(target, group)
        diff_results = {}
        if not run_diff:
            analysis_notes.append(
                "pairwise differential accessibility skipped (run_diff=false); "
                "marker_regions is cluster-vs-rest. Set run_diff=true for cluster-vs-cluster diff_test"
            )
        elif _labels is not None:
            _labels = _labels.astype(str)
            _levels = sorted_levels(_labels)
            for _g1, _g2 in combinations(_levels, 2):
                _idx1 = np.flatnonzero(_labels == _g1).tolist()
                _idx2 = np.flatnonzero(_labels == _g2).tolist()
                if not _idx1 or not _idx2:
                    continue
                try:
                    _df = snap.tl.diff_test(
                        target, cell_group1=_idx1, cell_group2=_idx2,
                        min_log_fc=diff_min_log_fc, min_pct=diff_min_pct, direction="both",
                    )
                    _key = f"{_g1}_vs_{_g2}"
                    diff_results[_key] = int(_df.height)
                    if _df.height:
                        _df.write_csv(os.path.join(tables_dir, f"diff_{_key}.csv"))
                except Exception as error:
                    info(f"diff_test {_g1} vs {_g2} skipped: " + str(error))
            analysis_steps.append("diff-accessibility")
            print(f"differential accessibility pairs: {diff_results}")
    except Exception as error:
        analysis_notes.append("markers skipped: " + str(error))
        info("markers skipped: " + str(error))