# JSON contract for sc_py_annotate.
notes = []
if not panel_path and not mapping_path:
    notes.append("no --panel and no --mapping: nothing was annotated")
if annotation_scores is not None:
    unassigned = [entry["cluster"] for entry in annotation_scores if entry["best"] is None]
    if unassigned:
        notes.append(
            f"{len(unassigned)} cluster(s) scored below --min-score {min_score}: "
            + ", ".join(unassigned)
        )
    if any(s["genes_missing"] > 0 for entry in annotation_scores for s in entry["scores"]):
        notes.append(
            "some panel genes are absent from this object; check the gene symbol "
            "naming (human MT- vs mouse mt-)"
        )

emit(
    {
        "input": tool_input,
        "output": tool_output,
        "group_by": group_by,
        "label_column": label_column if annotation_assignments is not None else None,
        "method": annotation_method,
        "panel": annotation_panel,
        "scores": annotation_scores,
        "assignments": annotation_assignments,
        "automated": annotation_automated,
        "figures": [f for f in [figure_paths.get("annotation-dotplot")] if f],
        "notes": notes,
    }
)