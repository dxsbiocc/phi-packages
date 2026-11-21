# JSON contract for scatac_markers.
_mk = markers or {}
total = sum(len(v) for v in _mk.values())
groups = [{"cluster": str(k), "n_markers": len(v)} for k, v in _mk.items()]

notes = list(analysis_notes)
notes.append(
    "clusters were defined from the same data (double dipping); treat marker rankings as evidence, not significance"
)
if total == 0:
    notes.append("no marker region passed the thresholds; relax marker_pvalue or check that clusters differ at all")
if not diff_results:
    notes.append("diff_accessibility is empty: pairwise diff_test is opt-in (run_diff=true); marker_regions is cluster-vs-rest")

emit(
    {
        "input": tool_input,
        "output": tool_output,
        "group_by": group_by,
        "marker_pvalue": marker_pvalue,
        "total_markers": total,
        "groups": groups,
        "tables": marker_tables,
        "diff_accessibility": diff_results,
        "diff_min_log_fc": diff_min_log_fc,
        "diff_min_pct": diff_min_pct,
        "notes": notes,
    }
)