# JSON contract for st_py_deconvolve. A skipped run carries {"skipped": reason}.
emit(
    {
        "input": tool_input,
        "output": tool_output,
        "deconvolution": deconv_results,
        "figures": [f for f in (figure_paths.get("deconv-proportions-spatial"),) if f],
        "notes": analysis_notes,
    }
)