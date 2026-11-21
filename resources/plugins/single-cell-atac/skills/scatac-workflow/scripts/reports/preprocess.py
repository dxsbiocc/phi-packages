# JSON contract for scatac_preprocess.
steps = [
    s
    for s in analysis_steps
    if s in ("tile-matrix", "select-features", "doublet-filter", "spectral", "knn", "clusters")
]

_sel = var_col(adata, "selected")
n_selected = int(np.asarray(_sel).sum()) if _sel is not None else adata.n_vars

cluster_summary = []
for column in cluster_columns:
    values = obs_col(adata, column)
    if values is None:
        continue
    names, counts = np.unique(values.astype(str), return_counts=True)
    order = sorted_levels(names)
    size_of = {str(n): int(c) for n, c in zip(names, counts)}
    cluster_summary.append(
        {
            "column": column,
            "clusters": int(len(names)),
            "sizes": [{"cluster": c, "cells": size_of[c]} for c in order],
        }
    )

emit(
    {
        "input": tool_input,
        "output": tool_output,
        "genome": genome_info,
        "steps": steps,
        "parameters": {
            "bin_size": bin_size,
            "n_features": n_features,
            "n_comps": n_comps,
            "n_neighbors": n_neighbors,
            "resolution": resolution,
            "doublet_filter": doublet_filter,
            "seed": seed,
        },
        "n_selected_features": n_selected,
        "clusters": cluster_summary,
        "embeddings": embeddings,
        "figures": [f for f in (figure_paths.get("dimplot-clusters"),) if f],
        "notes": analysis_notes,
    }
)