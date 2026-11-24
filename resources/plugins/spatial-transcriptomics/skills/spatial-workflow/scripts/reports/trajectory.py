# JSON contract for st_py_trajectory: PAGA + diffusion pseudotime mapped onto
# tissue coordinates, plus optional gene-program spatial scoring.
emit(
    {
        "input": tool_input,
        "output": tool_output,
        "group_by": group_by,
        "trajectory": trajectory_results,
        "gene_scores": score_results,
        "figures": [
            f
            for f in (
                figure_paths.get("trajectory-paga"),
                figure_paths.get("trajectory-spatial-dpt"),
                figure_paths.get("trajectory-gene-scores"),
            )
            if f
        ],
        "notes": analysis_notes,
    }
)
