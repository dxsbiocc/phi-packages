# Neighbourhood enrichment asks: do two domains sit next to each other more often
# than chance? It counts cluster-cluster edges in the spatial graph, then shuffles
# the labels n_perms times to build a null distribution; the z-score says how far
# the observed count is from that null. A high positive z-score on the diagonal
# means a domain is self-adjacent (a coherent region); off-diagonal highs mark
# interfaces between domains (e.g. tumour boundary next to immune infiltrate).
#
# Co-occurrence is the distance-resolved view: the probability of finding domain B
# within a growing radius of domain A, divided by B's overall frequency. A curve
# above 1 at short distances means A and B are locally mixed; the distance where it
# crosses 1 is the scale of their spatial association.
import squidpy as sq
import matplotlib.pyplot as plt

if "spatial_connectivities" not in adata.obsp:
    raise ValueError("no spatial graph; run the spatial-neighbors step first")
if group_by not in adata.obs.columns:
    raise ValueError(f"group column '{group_by}' not in obs")
if not isinstance(adata.obs[group_by].dtype, pd.CategoricalDtype):
    adata.obs[group_by] = adata.obs[group_by].astype("category")
categories = list(adata.obs[group_by].cat.categories)
if len(categories) < 2:
    info(f"only {len(categories)} group in {group_by}; neighbourhood enrichment needs >= 2, skipping")
    analysis_notes.append("nhood skipped: fewer than 2 groups")
else:
    sq.gr.nhood_enrichment(adata, cluster_key=group_by, n_perms=n_perms, n_jobs=1, show_progress_bar=False, seed=seed)
    sq.gr.co_occurrence(adata, cluster_key=group_by, spatial_key=coord_key, interval=cooc_interval)
    analysis_steps.append("nhood")

    zscore = adata.uns[f"{group_by}_nhood_enrichment"]["zscore"]
    pairs = []
    for i, a in enumerate(categories):
        for j, b in enumerate(categories):
            if i != j:
                pairs.append({"from": str(a), "to": str(b), "zscore": round(float(zscore[i, j]), 3)})
    pairs.sort(key=lambda p: -p["zscore"])
    nhood_results = {
        "cluster_key": group_by,
        "n_perms": n_perms,
        "groups": [str(c) for c in categories],
        "top_enriched_pairs": pairs[:10],
    }
    occ = adata.uns[f"{group_by}_co_occurrence"]
    cooc_results = {
        "cluster_key": group_by,
        "interval": cooc_interval,
        "distance_range": [float(occ["interval"][1]), float(occ["interval"][-1])],
        "shape": list(occ["occ"].shape),
    }
    print("top enriched domain pairs (z-score):")
    for p in pairs[:6]:
        print(f"  {p['from']} -> {p['to']}: {p['zscore']}")

    sq.pl.nhood_enrichment(adata, cluster_key=group_by, annotate=True, cmap="viridis")
    show_plot(plt.gcf(), "nhood-enrichment")
    # Co-occurrence panels grow with the number of domains; show a readable subset.
    show_clusters = [str(c) for c in categories[:4]]
    sq.pl.co_occurrence(adata, cluster_key=group_by, clusters=show_clusters)
    show_plot(plt.gcf(), "nhood-cooccurrence")