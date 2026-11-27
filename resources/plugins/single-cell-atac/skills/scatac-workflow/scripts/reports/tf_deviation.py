# JSON contract for scatac_tf_deviation.
_info = tfdev_info or {}
_tables = [p for p in (_info.get("output_deviations"), _info.get("output_variability")) if p]
emit(
    {
        "input": tool_input,
        "output": tool_output,
        **_info,
        "tables": _tables,
        "figures": [f for f in (figure_paths.get("tf-deviation-summary"),) if f],
        "notes": analysis_notes,
        "steps": list(analysis_steps),
    }
)
