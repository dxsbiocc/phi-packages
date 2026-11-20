# Map the expression-space clusters back onto tissue. This is the payoff of
# spatial data: a cluster is a spatial domain when its spots sit together in the
# tissue, and a scattered cluster is a cell state spread across regions. The map
# uses obsm[coord_key] directly (sc.pl.embedding), so it works for Visium spots,
# Stereo-seq bins and cellbin alike — no tissue image is required.
if coord_key not in adata.obsm:
    raise ValueError(f"coordinates '{coord_key}' not in obsm; cannot map to tissue")
if not cluster_columns:
    raise ValueError("no cluster column found; run the cluster step first")

column = cluster_columns[0]
show_plot(
    sc.pl.embedding(adata, basis=coord_key, color=column, legend_loc="on data", show=False),
    "spatial-clusters",
)

# A few top markers per domain, mapped in space, turn abstract cluster IDs into
# anatomy. Read these next to the cluster map: a marker that paints a known
# structure confirms the domain; one that paints nothing useful questions it.
if markers and markers.get("groups"):
    genes = []
    for entry in markers["groups"][:4]:
        genes.extend(entry["top_genes"][:2])
    genes = [g for g in dict.fromkeys(genes) if g in adata.var_names][:8]
    if genes:
        show_plot(
            sc.pl.embedding(
                adata,
                basis=coord_key,
                color=genes,
                layer="lognorm" if "lognorm" in adata.layers else None,
                show=False,
            ),
            "spatial-markers",
        )
analysis_steps.append("spatial-plots")