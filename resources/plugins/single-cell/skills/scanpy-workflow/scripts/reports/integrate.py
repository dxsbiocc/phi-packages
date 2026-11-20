# JSON contract for sc_py_integrate.
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

steps = [integration_info["method"], "neighbors", "clusters"]
if "umap" in embeddings:
    steps.append("umap")

emit(
    {
        "input": tool_input,
        "output": tool_output,
        "method": integration_info["method"],
        "batch_key": integration_info["batch_key"],
        "batches": integration_info["batches"],
        "embedding": integration_info["embedding"],
        "steps": steps,
        "clusters": cluster_summary,
        "figures": [
            f for f in (figure_paths.get("batch-before"), figure_paths.get("batch-after")) if f
        ],
        "notes": integration_info["notes"],
    }
)