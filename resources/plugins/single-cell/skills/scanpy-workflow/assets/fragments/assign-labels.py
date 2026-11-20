# Write the labels back. Edit this by hand if your reading of the dot plot
# disagrees with the automatic call — that judgement is the analysis. A reference
# panel forces a label onto every cluster, so a cluster with no evidence stays
# "unassigned" rather than being mislabelled.
import json

if mapping_path:
    if not os.path.exists(mapping_path):
        raise ValueError("mapping not found: " + mapping_path)
    try:
        with open(mapping_path) as handle:
            mapping = json.load(handle)
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError("cannot parse mapping JSON: " + str(error))
    mapping = {str(k).replace("-", "_"): str(v) for k, v in mapping.items()}
    present = sorted_levels(adata.obs[group_by])
    missing = [cluster for cluster in present if cluster not in mapping]
    if missing:
        raise ValueError(
            "mapping is missing clusters: " + ", ".join(missing)
            + "; label every cluster or omit --mapping"
        )
    adata.obs[label_column] = adata.obs[group_by].astype(str).map(mapping).to_numpy()
    info(f"wrote {label_column} from mapping")
    annotation_assignments = [{"cluster": c, "label": mapping[c]} for c in present]
elif annotation_scores is not None:
    labels = {
        entry["cluster"]: entry["best"]["cell_type"] if entry["best"] else "unassigned"
        for entry in annotation_scores
    }
    adata.obs[label_column] = adata.obs[group_by].astype(str).map(labels).to_numpy()
    annotation_assignments = [
        {
            "cluster": entry["cluster"],
            "label": entry["best"]["cell_type"] if entry["best"] else "unassigned",
            "score": entry["best"]["score"] if entry["best"] else None,
        }
        for entry in annotation_scores
    ]
    info(f"wrote {label_column} from the highest-scoring panel per cluster")
    if (adata.obs[label_column] == "unassigned").all():
        raise ValueError("No cluster could be scored against the panel; fix the identifiers before labelling.")

if annotation_assignments is not None:
    print(adata.obs[label_column].value_counts())
    if "X_umap" in adata.obsm:
        show_plot(sc.pl.umap(adata, color=label_column, legend_loc="on data", show=False), "dimplot-celltypes")