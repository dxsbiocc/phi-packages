# JSON contract for scatac_annotate.
emit(
    {
        "input": tool_input,
        "output": tool_output,
        "group_by": group_by,
        "label_column": label_column,
        "method": annotation_automated,
        "panel": annotation_panel,
        "assignments": annotation_assignments,
        "scores": annotation_scores,
        "gene_activity_present": gene_adata is not None,
        "notes": analysis_notes
        + (
            []
            if annotation_automated
            else ["no annotation produced; pass gene_anno for panel scoring, a reference_path, or a mapping_path"]
        ),
    }
)