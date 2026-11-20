# Edit this panel to your tissue; it is the evidence the labels rest on. A label
# needs positive markers AND the absence of markers of competing types, so read
# the dot plot below before accepting the automatic call. The built-in panel is a
# generic example — spatial regions often need tissue-specific markers, so pass
# panel_path with a curated JSON (cell type -> array of gene symbols).
if group_by not in adata.obs.columns:
    raise ValueError(
        f"group column '{group_by}' not in obs; available: " + ", ".join(adata.obs.columns)
    )
if panel_path:
    import json

    if not os.path.exists(panel_path):
        raise ValueError("panel not found: " + panel_path)
    try:
        with open(panel_path) as handle:
            parsed = json.load(handle)
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError("cannot parse panel JSON: " + str(error))
    if not isinstance(parsed, dict):
        raise ValueError("panel JSON must be an object mapping cell type -> array of gene symbols")
    panel = {}
    for cell_type, genes in parsed.items():
        if not isinstance(genes, list) or not all(isinstance(g, str) for g in genes):
            raise ValueError("each panel entry must be an array of gene symbols")
        panel[cell_type] = genes
else:
    panel = {
        "T_cell": ["CD3D", "CD3E", "TRAC", "IL7R", "LTB"],
        "B_cell": ["MS4A1", "CD79A", "CD79B", "CD19"],
        "Plasma": ["MZB1", "XBP1", "JCHAIN", "SDC1"],
        "NK": ["NKG7", "KLRD1", "GNLY", "PRF1"],
        "Myeloid": ["LYZ", "S100A8", "S100A9", "CTSS"],
        "Dendritic": ["FCER1A", "CLEC10A", "CD1C", "ITGAX"],
        "Fibroblast": ["COL3A1", "PDGFRA", "DCN", "LAMA2"],
        "Epithelial": ["EPCAM", "KRT8", "KRT18", "CDH1"],
        "Endothelial": ["PECAM1", "VWF", "CDH5", "CLDN5"],
    }

expression = adata.layers["lognorm"] if "lognorm" in adata.layers else adata.X
all_genes = list(adata.var_names)
gene_index = {gene: i for i, gene in enumerate(all_genes)}
missing = sorted({g for genes in panel.values() for g in genes} - set(gene_index))
print(f"panel genes missing from this object: {len(missing)}")
if missing:
    print(missing)
usable = {ct: [g for g in genes if g in gene_index] for ct, genes in panel.items()}
usable = {ct: genes for ct, genes in usable.items() if genes}
if not usable:
    raise ValueError(
        "None of the panel genes occur in this object. Check the identifier type: "
        f"var_names hold {all_genes[:3]}... Human symbols use MT-, mouse mt-; "
        "Ensembl IDs need a different panel."
    )

# Score = panel mean minus the mean of everything else, so panel sizes stay
# comparable. This is a heuristic ranking, not a trained classifier: a spatial
# spot mixes several cells, so a high score means "this region is enriched for",
# not "this spot is one cell of that type".
annotation_panel = list(panel.keys())
annotation_scores = []
for group in sorted_levels(adata.obs[group_by]):
    mask = (adata.obs[group_by].astype(str) == group).to_numpy()
    units = int(mask.sum())
    if units == 0:
        continue
    mean = np.asarray(expression[mask].mean(axis=0)).ravel()
    total_mean = float(mean.sum())
    scores = []
    for cell_type, genes in panel.items():
        found = [gene_index[g] for g in genes if g in gene_index]
        if not found:
            scores.append(
                {"cell_type": cell_type, "score": None, "genes_found": 0, "genes_missing": len(genes)}
            )
            continue
        panel_mean = float(mean[found].mean())
        rest_count = len(all_genes) - len(found)
        rest_mean = (total_mean - float(mean[found].sum())) / rest_count if rest_count else 0.0
        scores.append(
            {
                "cell_type": cell_type,
                "score": round(panel_mean - rest_mean, 4),
                "panel_mean": round(panel_mean, 4),
                "genes_found": len(found),
                "genes_missing": len(genes) - len(found),
            }
        )
    scored = sorted((s for s in scores if s["score"] is not None), key=lambda s: -s["score"])
    best = scored[0] if scored and scored[0]["score"] >= min_score else None
    annotation_scores.append({"cluster": group, "units": units, "best": best, "scores": scores})

best_table = pd.DataFrame(
    [
        {
            "cluster": entry["cluster"],
            "label": entry["best"]["cell_type"] if entry["best"] else "unassigned",
            "score": entry["best"]["score"] if entry["best"] else None,
        }
        for entry in annotation_scores
    ]
)
print(best_table.to_string(index=False))

plot_genes = list(dict.fromkeys(g for genes in usable.values() for g in genes))[:60]
if plot_genes:
    sc.tl.dendrogram(adata, groupby=group_by)
    show_plot(
        sc.pl.dotplot(
            adata,
            var_names=plot_genes,
            groupby=group_by,
            layer="lognorm" if "lognorm" in adata.layers else None,
            standard_scale="var",
            show=False,
        ),
        "annotation-dotplot",
    )