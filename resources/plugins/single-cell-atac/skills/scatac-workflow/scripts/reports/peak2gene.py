# JSON contract for scatac_peak2gene.
_info = peak2gene_info or {}
_tables = [p for p in (_info.get("output_links"),) if p]
emit(
    {
        "input": tool_input,
        "tables": _tables,
        **_info,
        "figures": [f for f in (figure_paths.get("peak2gene-summary"),) if f],
        "notes": analysis_notes,
        "steps": list(analysis_steps),
    }
)
