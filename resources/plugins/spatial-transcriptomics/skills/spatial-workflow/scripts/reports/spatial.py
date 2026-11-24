# JSON contract for st_py_spatial: SVG, neighbourhood enrichment, co-occurrence,
# interaction matrix, centrality scores, Ripley's L/F/G and cross-sample
# differential spatial patterns. A skipped section carries a {"skipped": reason}
# object so its absence is never read as a clean result.
emit(
    {
        "input": tool_input,
        "output": tool_output,
        "group_by": group_by,
        "svg": svg_results,
        "nhood_enrichment": nhood_results,
        "co_occurrence": cooc_results,
        "interaction_matrix": interaction_results,
        "centrality": centrality_results,
        "ripley": ripley_results,
        "spatial_diff": spatial_diff_results,
        "figures": [
            f
            for f in (
                figure_paths.get("svg-top"),
                figure_paths.get("svg-top-gene-spatial"),
                figure_paths.get("nhood-enrichment"),
                figure_paths.get("nhood-cooccurrence"),
                figure_paths.get("nhood-centrality"),
                figure_paths.get("nhood-ripley"),
                figure_paths.get("spatial-diff-moran"),
            )
            if f
        ],
        "notes": analysis_notes,
    }
)
