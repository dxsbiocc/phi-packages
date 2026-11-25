# JSON contract for sc_py_communication.
emit(
    {
        "input": tool_input,
        "output": tool_output,
        **(ccc_info or {}),
        "figures": [
            f
            for f in (
                figure_paths.get("communication-heatmap"),
                figure_paths.get("communication-dotplot"),
            )
            if f
        ],
        "notes": analysis_notes
        + [
            "mRNA ligand-receptor co-expression indicates signaling potential, not physical cell-cell contact",
            "pass --interactions_path with an OmniPath/CellPhoneDB TSV for exhaustive receptorome screening, or use liana in a project environment for multi-method consensus",
        ],
    }
)
