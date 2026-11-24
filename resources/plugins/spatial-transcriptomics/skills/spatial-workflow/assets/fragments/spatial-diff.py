# Cross-sample differential spatial patterns ask: does a gene's spatial structure
# differ between conditions (e.g. disease vs control, region A vs B)? A gene can be
# strongly spatially organised in one condition and diffuse in another — that change
# is the spatial phenotype, invisible to ordinary differential expression which
# ignores location.
#
# For each condition we subset the spots, rebuild their OWN spatial graph (their
# layout differs), and recompute Moran's I. Comparing I across conditions ranks
# genes by how much their spatial organisation changes. We also report each
# condition's cluster composition, so a shift in domain proportions can be read
# alongside the shift in spatial structure. This is a descriptive screen, not a
# significance test of the difference itself.
import squidpy as sq
import matplotlib.pyplot as plt
import seaborn as sns

if not group_column:
    info("no group_column; skipping cross-sample differential. Set it to a condition/sample obs column.")
    spatial_diff_results = {"skipped": "no group_column provided"}
    analysis_notes.append("spatial-diff skipped: no group_column")
elif group_column not in adata.obs.columns:
    raise ValueError(f"group column '{group_column}' not in obs; available: " + ", ".join(adata.obs.columns))
else:
    conditions = sorted_levels(adata.obs[group_column].astype(str))
    layer = "lognorm" if "lognorm" in adata.layers else None
    moran_by_condition = {}
    sig_by_condition = {}
    composition = {}
    used = []
    for condition in conditions:
        mask = (adata.obs[group_column].astype(str) == condition).to_numpy()
        sub = adata[mask].copy()
        if sub.n_obs < min_units:
            info(f"condition '{condition}' has {sub.n_obs} units (< min_units={min_units}); skipping")
            continue
        if hasattr(sq.gr, "spatial_neighbors_knn"):
            sq.gr.spatial_neighbors_knn(sub, spatial_key=coord_key, n_neighs=n_neighs)
        else:
            sq.gr.spatial_neighbors(sub, spatial_key=coord_key, coord_type="generic", n_neighs=n_neighs)
        sq.gr.spatial_autocorr(
            sub, mode="moran", genes=None, layer=layer,
            corr_method="fdr_bh", n_perms=None, use_raw=False, seed=seed,
        )
        frame = sub.uns["moranI"]
        moran_by_condition[condition] = frame["I"]
        sig_by_condition[condition] = set(frame.index[frame["pval_norm_fdr_bh"] < 0.05])
        used.append(condition)
        if group_by in sub.obs.columns:
            counts = sub.obs[group_by].astype(str).value_counts(normalize=True)
            composition[condition] = {str(k): round(float(v), 4) for k, v in counts.items()}

    if len(used) < 2:
        info("fewer than 2 conditions with enough units; skipping cross-sample differential")
        spatial_diff_results = {"skipped": "fewer than 2 usable conditions"}
        analysis_notes.append("spatial-diff skipped: fewer than 2 usable conditions")
    else:
        analysis_steps.append("spatial-diff")
        table = pd.DataFrame(moran_by_condition)
        table["range"] = table.max(axis=1) - table.min(axis=1)
        table = table.sort_values("range", ascending=False)
        top_genes = list(table.index[:20])

        # Condition-specific SVGs: significant in exactly one condition.
        all_sig = set().union(*sig_by_condition.values())
        specific = {}
        for gene in all_sig:
            hits = [c for c in used if gene in sig_by_condition[c]]
            if len(hits) == 1:
                specific.setdefault(hits[0], []).append(str(gene))
        spatial_diff_results = {
            "group_column": group_column,
            "conditions": used,
            "n_genes_compared": int(len(table)),
            "top_differential_genes": [
                {"gene": str(g), "moran_I": {c: round(float(table.loc[g, c]), 4) for c in used},
                 "range": round(float(table.loc[g, "range"]), 4)}
                for g in top_genes[:10]
            ],
            "condition_specific_svg": {c: genes[:20] for c, genes in specific.items()},
            "composition": composition,
        }
        print(f"compared Moran's I across {len(used)} conditions for {len(table)} genes")
        print("top genes by spatial-structure change:")
        print(table[used].head(10).round(3).to_string())

        plot_table = table.loc[top_genes, used]
        fig, ax = plt.subplots(figsize=(max(4.0, len(used) * 1.3), max(4.0, len(top_genes) * 0.3)))
        sns.heatmap(plot_table, annot=True, fmt=".2f", cmap="viridis", ax=ax, cbar_kws={"label": "Moran's I"})
        ax.set_title("Spatial structure (Moran's I) by condition")
        show_plot(fig, "spatial-diff-moran")