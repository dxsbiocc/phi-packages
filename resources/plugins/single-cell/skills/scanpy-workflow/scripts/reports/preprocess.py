# JSON contract for sc_py_preprocess.
steps = [
    s
    for s in analysis_steps
    if s in ("normalize", "variable-features", "regress", "scale", "pca", "neighbors", "clusters")
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
                {"cluster": name, "cells": int(counts_table[name])}
                for name in sorted_levels(counts_table.index)
            ],
        }
    )

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
            "seed": seed,
        },
        "top_variable_features": top_variable,
        "variance": pca_variance,
        "clusters": cluster_summary,
        "embeddings": embeddings,
        "figures": [
            f for f in (figure_paths.get("dimplot-clusters"), figure_paths.get("pca-elbow")) if f
        ],
        "notes": analysis_notes,
    }
)