# JSON contract for scatac_peaks.
emit(
    {
        "input": tool_input,
        "output": tool_output,
        "group_by": group_by,
        "peak_qvalue": peak_qvalue,
        "peak_half_width": peak_half_width,
        "call_broad_peaks": call_broad_peaks,
        "peaks_per_cluster": peak_tables,
        "merged_peaks": merged_peaks,
        "steps": [s for s in analysis_steps if s in ("peaks",)],
        "notes": analysis_notes
        + (
            []
            if merged_peaks
            else ["no merged peak set; peak calling needs fragments and a cluster column, and a permissive enough qvalue"]
        ),
    }
)