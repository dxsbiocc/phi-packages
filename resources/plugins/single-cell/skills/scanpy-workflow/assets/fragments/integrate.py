# Batch effects are technical differences between samples (donor, chemistry,
# run). Integration aligns the shared biology so clusters are not one-per-sample.
# Harmony corrects the PCA embedding; ComBat corrects log-normalized expression
# and needs a fresh PCA. Judge the result on the before/after UMAPs: samples
# should mix within each cell type, and marker genes must survive — over-
# correction erases real biology (e.g. patient-specific tumour states).
if batch_key not in adata.obs.columns:
    print(f"No '{batch_key}' column in this object; skipping integration.")
elif "X_pca" not in adata.obsm:
    raise ValueError("no X_pca embedding; run sc_py_preprocess before sc_py_integrate")
else:
    batches = [
        {"batch": str(name), "cells": int(count)}
        for name, count in adata.obs[batch_key].value_counts().items()
    ]
    notes = []
    if batch_method == "harmony":
        # scanpy 1.13 resolves batch as an accessor string, not a bare column name.
        sc.pp.harmony_integrate(
            adata, f"obs.{batch_key}", use="obsm.X_pca", out="obsm.X_pca_harmony", rng=seed
        )
        adata.obsm["X_pca_integrated"] = adata.obsm["X_pca_harmony"].copy()
        print(f"harmony corrected the PCA axis on '{batch_key}'")
    elif batch_method == "combat":
        # combat models log-normalized expression, not the scaled matrix.
        if "lognorm" in adata.layers:
            source = adata.copy()
            source.X = adata.layers["lognorm"].copy()
        else:
            source = adata
        sc.pp.combat(source, batch_key)
        sc.tl.pca(
            source,
            n_comps=n_pcs,
            mask="var.highly_variable" if "highly_variable" in source.var.columns else None,
            svd_solver="arpack",
            rng=seed,
        )
        adata.obsm["X_pca_integrated"] = source.obsm["X_pca"].copy()
        notes.append("combat ran on the pre-scale lognorm expression snapshot")
        print(f"combat corrected expression on '{batch_key}'; PCA recomputed")
    elif batch_method == "bbknn":
        if not optional_module("bbknn"):
            raise ValueError(
                "bbknn is not installed in plugin:scanpy (no python 3.13 build on bioconda); "
                "use --method harmony or --method combat, or install bbknn into a project "
                "environment with python <= 3.12"
            )
        sc.external.pp.bbknn(adata, batch_key=batch_key)
        adata.obsm["X_pca_integrated"] = adata.obsm["X_pca"].copy()
        notes.append("bbknn corrects the neighbour graph, not the embedding")
    else:
        raise ValueError("--method must be harmony, combat or bbknn, got " + batch_method)
    analysis_steps.append(batch_method)
    analysis_notes.extend(notes)
    integration_info = {
        "method": batch_method,
        "batch_key": batch_key,
        "batches": batches,
        "embedding": "X_pca_integrated",
        "notes": notes,
    }