# JSON contract for sc_py_markers.
if markers["total"] == 0:
    emit(
        {
            "input": tool_input,
            "output": output_object,
            "group_by": group_by,
            "test": test_use,
            "total_markers": 0,
            "groups": [],
            "notes": [
                "no marker passed the thresholds; lower --logfc-threshold or --min-pct, "
                "or check that clusters differ at all"
            ],
        }
    )
    sys.exit(0)

notes = []
if only_positive:
    notes.append("only positive markers are reported (--only-positive)")
if markers["snapshot"]:
    notes.append("markers tested on the pre-scale lognorm snapshot")
else:
    notes.append("no lognorm snapshot; markers tested on the current X")

emit(
    {
        "input": tool_input,
        "output": tool_output,
        "group_by": group_by,
        "test": test_use,
        "total_markers": markers["total"],
        "groups": markers["groups"],
        "tables": marker_tables,
        "figures": [
            f
            for f in [figure_paths.get("markers-heatmap-" + re.sub(r"[^A-Za-z0-9._-]", "_", group_by))]
            if f
        ],
        "notes": notes,
    }
)