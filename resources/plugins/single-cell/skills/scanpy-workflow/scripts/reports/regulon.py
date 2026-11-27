# JSON contract for sc_py_regulon.
emit(
    {
        "input": tool_input,
        "output": tool_output,
        **(reg_info or {}),
        "figures": [f for f in (figure_paths.get("regulon-heatmap"),) if f],
        "notes": analysis_notes
        + [
            "expression-based regulon activity reflects target co-expression, not direct chromatin occupancy",
            "pass --regulons_path (CollecTRI/DoRothEA TSV or .gmt) to test an expanded TF panel",
        ],
    }
)
