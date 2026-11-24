# Neighbourhood enrichment asks: do two domains sit next to each other more often
# than chance? It counts cluster-cluster edges in the spatial graph, then shuffles
# the labels n_perms times to build a null distribution; the z-score says how far
# the observed count is from that null. A high positive z-score on the diagonal
# means a domain is self-adjacent (a coherent region); off-diagonal highs mark
# interfaces between domains (e.g. tumour boundary next to immune infiltrate).
#
# Co-occurrence is the distance-resolved view: the probability of finding domain B
# within a growing radius of domain A, divided by B's overall frequency.
#
# Three complementary topology views run on the same spatial graph:
# - Interaction matrix: row-normalized contact frequencies between domains.
# - Centrality scores: degree, clustering, and closeness centrality of each domain
#   on the tissue graph (distinguishing core vs dispersed/boundary domains).
# - Ripley's L/F/G: point-pattern clustering across spatial scales; L(r) - r > 0
#   indicates clustering at distance r, < 0 indicates dispersion.
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
    sq.gr.co_occurrence(adata, cluster_key=group_by, spatial_key=coord_key, interval=cooc_interval, n_jobs=1, show_progress_bar=False)
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

    # Interaction matrix: normalized domain-to-domain contact frequencies.
    try:
        imat = sq.gr.interaction_matrix(adata, cluster_key=group_by, normalized=True, copy=True)
        adata.uns[f"{group_by}_interactions"] = imat
        contact_pairs = []
        for i, a in enumerate(categories):
            for j, b in enumerate(categories):
                if i != j:
                    contact_pairs.append(
                        {"from": str(a), "to": str(b), "proportion": round(float(imat[i, j]), 4)}
                    )
        contact_pairs.sort(key=lambda p: -p["proportion"])
        interaction_results = {
            "cluster_key": group_by,
            "normalized": True,
            "top_contact_pairs": contact_pairs[:10],
        }
    except Exception as error:
        info(f"interaction_matrix skipped: {error}")
        interaction_results = {"skipped": str(error)}
        analysis_notes.append(f"interaction_matrix skipped: {error}")

    # Spatial graph centrality per domain.
    try:
        sq.gr.centrality_scores(adata, cluster_key=group_by, n_jobs=1, show_progress_bar=False)
        cframe = adata.uns[f"{group_by}_centrality_scores"]
        centrality_results = {
            "cluster_key": group_by,
            "scores": {
                str(idx): {
                    str(col): round(float(val), 4)
                    for col, val in row.items()
                    if pd.notna(val)
                }
                for idx, row in cframe.iterrows()
            },
        }
        sq.pl.centrality_scores(adata, cluster_key=group_by)
        show_plot(plt.gcf(), "nhood-centrality")
    except Exception as error:
        info(f"centrality_scores skipped: {error}")
        centrality_results = {"skipped": str(error)}
        analysis_notes.append(f"centrality_scores skipped: {error}")

    # Ripley's point-pattern statistic (L, F, or G).
    try:
        _ripley_mode = ripley_mode if ripley_mode in ("L", "F", "G") else "L"
        sq.gr.ripley(
            adata,
            cluster_key=group_by,
            mode=_ripley_mode,
            spatial_key=coord_key,
            n_simulations=int(ripley_simulations),
            seed=seed,
        )
        rkey = f"{group_by}_ripley_{_ripley_mode}"
        rdata = adata.uns[rkey]
        stat_df = rdata[f"{_ripley_mode}_stat"]
        per_cluster = {}
        for cluster_name, grp in stat_df.groupby(group_by, observed=True):
            bins = grp["bins"].to_numpy(dtype=float)
            stats = grp["stats"].to_numpy(dtype=float)
            if _ripley_mode == "L":
                excess = stats - bins
                per_cluster[str(cluster_name)] = {
                    "max_excess_L_minus_r": round(float(np.nanmax(excess)), 4),
                    "min_excess_L_minus_r": round(float(np.nanmin(excess)), 4),
                    "pattern": "clustered" if float(np.nanmean(excess)) > 0 else "dispersed",
                }
            else:
                per_cluster[str(cluster_name)] = {
                    "max_stat": round(float(np.nanmax(stats)), 4),
                    "mean_stat": round(float(np.nanmean(stats)), 4),
                }
        ripley_results = {
            "cluster_key": group_by,
            "mode": _ripley_mode,
            "n_simulations": int(ripley_simulations),
            "per_cluster": per_cluster,
        }
        sq.pl.ripley(adata, cluster_key=group_by, mode=_ripley_mode)
        show_plot(plt.gcf(), "nhood-ripley")
    except Exception as error:
        info(f"ripley skipped: {error}")
        ripley_results = {"skipped": str(error)}
        analysis_notes.append(f"ripley skipped: {error}")
