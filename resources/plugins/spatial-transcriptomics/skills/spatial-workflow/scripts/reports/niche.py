# JSON contract for st_py_niche: spatial niche segmentation and optional
# cross-condition/stage niche dynamics.
emit(
    {
        "input": tool_input,
        "output": tool_output,
        "group_by": group_by,
        "niche": niche_results,
        "figures": [
            f
            for f in (
                figure_paths.get("niche-spatial"),
                figure_paths.get("niche-composition"),
                figure_paths.get("niche-condition-dynamics"),
            )
            if f
        ],
        "notes": analysis_notes,
    }
)
