# A spatially variable gene (SVG) has expression that depends on where a spot is,
# not just on which cluster it belongs to. Moran's I measures global spatial
# autocorrelation: I near +1 means neighbouring spots have similar expression
# (a smooth spatial pattern), I near 0 means no spatial structure, I < 0 means
# neighbours differ (checkerboard). Geary's C is the complementary ratio: C near 0
# is strong positive autocorrelation, C near 1 is none. Both read the spatial
# adjacency graph, so spatial-neighbors must run first.
#
# genes=None lets squidpy use the highly-variable genes when present (fast, and
# where spatial signal usually lives), otherwise all genes. p-values are analytic
# (normality assumption) — no permutations — which is the standard SVG screen;
# corrected across genes with a statsmodels method (fdr_bh by default).
import squidpy as sq
import matplotlib.pyplot as plt

if "spatial_connectivities" not in adata.obsp:
    raise ValueError("no spatial graph; run the spatial-neighbors step before SVG")
if svg_method not in ("moran", "geary"):
    raise ValueError("svg_method must be moran or geary")
corr = {"benjamini-hochberg": "fdr_bh", "bh": "fdr_bh"}.get(svg_corr_method, svg_corr_method)
layer = "lognorm" if "lognorm" in adata.layers else None

sq.gr.spatial_autocorr(
    adata,
    mode=svg_method,
    genes=None,
    layer=layer,
    corr_method=corr,
    n_perms=None,
    use_raw=False,
    seed=seed,
)
analysis_steps.append("svg")

key = "moranI" if svg_method == "moran" else "gearyC"
stat = "I" if svg_method == "moran" else "C"
frame = adata.uns[key]
pcol = f"pval_norm_{corr}"
significant = frame[frame[pcol] < 0.05] if pcol in frame else frame
# squidpy already sorts: Moran descending (high I first), Geary ascending (low C first).
top = frame.head(20)
svg_results = {
    "method": svg_method,
    "statistic": stat,
    "corr_method": corr,
    "n_genes": int(len(frame)),
    "n_significant": int(len(significant)),
    "top": [
        {"gene": str(g), stat: round(float(row[stat]), 4), "pval_adj": round(float(row[pcol]), 6)}
        for g, row in top.iterrows()
    ],
}
print(f"{svg_method}: {len(significant)} of {len(frame)} genes spatially variable (adj p < 0.05)")
print(top[[stat, pcol]].head(10).to_string())

fig, ax = plt.subplots(figsize=(6, max(3.0, len(top) * 0.28)))
ax.barh(np.arange(len(top)), top[stat].to_numpy(), color="#3b6fb6")
ax.set_yticks(np.arange(len(top)))
ax.set_yticklabels(top.index)
ax.invert_yaxis()
ax.set_xlabel(f"{svg_method} {stat}")
ax.set_title("Top spatially variable genes")
show_plot(fig, "svg-top")

if len(top) and coord_key in adata.obsm:
    show_plot(
        sc.pl.embedding(adata, basis=coord_key, color=[str(top.index[0])], layer=layer, show=False),
        "svg-top-gene-spatial",
    )