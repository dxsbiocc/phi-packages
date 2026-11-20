# A cluster (a spatial domain) is only a hypothesis until its markers say what it
# is. Wilcoxon ranks genes by how specifically they distinguish one group from the
# rest; the test runs on the pre-scale lognorm snapshot, because scaled values are
# z-scores and their "fold changes" are meaningless. p-values here are descriptive
# only: the clusters were defined from the same data (double dipping), so treat the
# ranking as evidence, not as significance.
import re

if group_by not in adata.obs.columns:
    raise ValueError(
        f"group column '{group_by}' not in obs; available: " + ", ".join(adata.obs.columns)
    )
if test_use not in ("wilcoxon", "t-test", "t-test_overestim_var", "logreg"):
    raise ValueError("test_use must be wilcoxon, t-test, t-test_overestim_var or logreg")

snapshot = "lognorm" in adata.layers
info(f"rank_genes_groups on {group_by} (test: {test_use}, lognorm={snapshot})")
sc.tl.rank_genes_groups(
    adata,
    groupby=group_by,
    method=test_use,
    layer="lognorm" if snapshot else None,
    n_genes=adata.n_vars,
    pts=True,
)

result = adata.uns["rank_genes_groups"]
groups = sorted_levels(result["names"].dtype.names)

per_group = []
marker_tables = []
total = 0
for group in groups:
    frame = pd.DataFrame(
        {
            "gene": result["names"][group],
            "score": result["scores"][group],
            "avg_log2FC": result["logfoldchanges"][group],
            "p_val": result["pvals"][group],
            "p_val_adj": result["pvals_adj"][group],
            "pts": result["pts"][group] if "pts" in result else np.nan,
        }
    )
    if only_positive:
        frame = frame[frame["avg_log2FC"] >= logfc_threshold]
    else:
        frame = frame[frame["avg_log2FC"].abs() >= logfc_threshold]
    if min_pct is not None and "pts" in result:
        frame = frame[frame["pts"] >= min_pct]
    frame = frame.sort_values(["p_val_adj", "avg_log2FC"], ascending=[True, False])
    frame.insert(1, "cluster", group)
    total += len(frame)
    if tables_dir and len(frame):
        os.makedirs(tables_dir, exist_ok=True)
        path = os.path.join(
            tables_dir, "markers-" + re.sub(r"[^A-Za-z0-9._-]", "_", group) + ".csv"
        )
        frame.to_csv(path, index=False)
        info("table " + path)
        marker_tables.append(path)
    top = frame.head(top_n)
    per_group.append(
        {
            "cluster": group,
            "markers": int(len(frame)),
            "top_genes": list(top["gene"]),
            "top_log2fc": [
                {"gene": row.gene, "avg_log2FC": round(float(row.avg_log2FC), 4)}
                for row in top.itertuples()
            ],
        }
    )
    print(f"cluster {group}: top {', '.join(top['gene'].head(5))}")

markers = {"total": int(total), "groups": per_group, "snapshot": snapshot}
show_plot(sc.pl.rank_genes_groups(adata, n_genes=top_n, sharey=False, show=False), "markers-rank")

if total:
    sc.tl.dendrogram(adata, groupby=group_by)
    show_plot(
        sc.pl.rank_genes_groups_heatmap(
            adata,
            n_genes=5,
            groupby=group_by,
            layer="lognorm" if snapshot else None,
            standard_scale="var",
            show=False,
        ),
        "markers-heatmap-" + re.sub(r"[^A-Za-z0-9._-]", "_", group_by),
    )