# JSON contract for st_py_preprocess.
steps = [
    s
    for s in analysis_steps
    if s in ("normalize", "variable-features", "scale", "pca", "neighbors", "clusters",
             "spatial-neighbors", "spatial-plots")
]

if "highly_variable_rank" in adata.var.columns:
    ranked = adata.var["highly_variable_rank"][adata.var["highly_variable"]]
    top_variable = list(ranked.sort_values().index[:10])
else:
    top_variable = list(adata.var_names[adata.var["highly_variable"]][:10])

cluster_summary = []
for column in cluster_columns:
    counts_table = adata.obs[column].value_counts()
    cluster_summary.append(
        {
            "column": column,
            "clusters": int(len(counts_table)),
            "sizes": [
                {"cluster": name, "units": int(counts_table[name])}
                for name in sorted_levels(counts_table.index)
            ],
        }
    )

spatial_graph = None
if "spatial_connectivities" in adata.obsp:
    graph = adata.obsp["spatial_connectivities"]
    degree = np.asarray(graph.sum(axis=1)).ravel()
    spatial_graph = {
        "method": neighbor_method,
        "edges": int(graph.nnz),
        "mean_degree": round(float(degree.mean()), 3),
    }

emit(
    {
        "input": tool_input,
        "output": tool_output,
        "steps": steps,
        "parameters": {
            "normalization": normalization,
            "scale_factor": scale_factor,
            "hvg_flavor": hvg_flavor,
            "n_variable_features": int(adata.var["highly_variable"].sum()),
            "n_pcs": pca_variance["pcs"],
            "dims_use": dims,
            "resolution": resolution,
            "neighbor_method": neighbor_method,
            "n_neighs": n_neighs,
            "seed": seed,
        },
        "top_variable_features": top_variable,
        "variance": pca_variance,
        "clusters": cluster_summary,
        "embeddings": embeddings,
        "spatial_graph": spatial_graph,
        "figures": [
            f
            for f in (
                figure_paths.get("dimplot-clusters"),
                figure_paths.get("pca-elbow"),
                figure_paths.get("hvg-mean-variance"),
                figure_paths.get("spatial-clusters"),
            )
            if f
        ],
        "notes": analysis_notes,
    }
)