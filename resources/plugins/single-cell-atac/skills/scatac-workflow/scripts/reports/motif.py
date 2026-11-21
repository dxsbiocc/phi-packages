# JSON contract for scatac_motif.
emit(
    {
        "input": tool_input,
        "output": tool_output,
        "motif_method": motif_method,
        "motif_pvalue": motif_pvalue,
        "results": motif_results,
        "steps": [s for s in analysis_steps if s in ("motif",)],
        "notes": analysis_notes
        + (
            []
            if motif_results
            else ["no motif enrichment produced; needs marker regions, a genome_fasta and a motif database (motif_db or cis_bp)"]
        ),
    }
)