# JSON contract for st_py_communication. A skipped run carries {"skipped": reason}.
emit(
    {
        "input": tool_input,
        "output": tool_output,
        "group_by": group_by,
        "communication": ccc_results,
        "figures": [f for f in (figure_paths.get("communication-ligrec"),) if f],
        "notes": analysis_notes,
    }
)