# JSON contract for scatac_integrate.
emit(
    {
        "input": tool_input,
        "output": tool_output,
        "integration": integration_info,
        "embeddings": embeddings,
        "steps": [s for s in analysis_steps if s in ("integrate",)],
        "notes": analysis_notes
        + (
            []
            if integration_info
            else ["no integration produced; needs a batch_column with at least two samples and a computed X_spectral"]
        ),
    }
)