# Annotation maps each cluster to a cell type. Three routes, in priority order:
#   1. mapping_path  — an explicit cluster -> label JSON the agent/curator decided.
#   2. panel scoring — score a marker panel on the GENE ACTIVITY matrix (gene_adata).
#      Gene activity is an accessibility proxy for expression, so a high score means
#      "this cluster's open chromatin is enriched for these genes' loci", not measured
#      RNA. It is a heuristic ranking, not a trained classifier — read the scores.
#   3. reference transfer — a scRNA reference (reference_path) co-embedded with gene
#      activity in a shared gene space, then kNN label transfer per cluster.
# All routes need either gene_anno (for gene activity) or an explicit mapping; without
# them the cell reports the gap. The chosen labels are written to obs[label_column].
import json

group = group_by if obs_has(adata, group_by) else (cluster_columns[0] if cluster_columns else None)
if group is None:
    raise ValueError("annotation needs a cluster column; run the cluster step first")

labels_by_cluster = {}
annotation_automated = None


def _load_panel():
    if panel_path:
        if not os.path.exists(panel_path):
            raise ValueError("panel not found: " + panel_path)
        with open(panel_path, encoding="utf-8") as handle:
            parsed = json.load(handle)
        if not isinstance(parsed, dict):
            raise ValueError("panel JSON must map cell type -> array of gene symbols")
        return {str(ct): [str(g) for g in genes] for ct, genes in parsed.items()}
    # Generic example panel; pass panel_path with tissue-specific markers for real data.
    return {
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


def _gene_memory():
    g = gene_adata.to_memory()
    if "gene_name" in g.var:
        g.var_names = pd.Index(g.var["gene_name"].astype(str))
        g.var_names_make_unique()
    return g


# 1) explicit mapping wins
if mapping_path:
    if not os.path.exists(mapping_path):
        raise ValueError("mapping not found: " + mapping_path)
    with open(mapping_path, encoding="utf-8") as handle:
        parsed = json.load(handle)
    labels_by_cluster = {str(k): str(v) for k, v in parsed.items()}
    annotation_automated = "mapping"

# 2) panel scoring on gene activity
if not labels_by_cluster and gene_adata is not None:
    g = _gene_memory()
    sc.pp.normalize_total(g, target_sum=1e4)
    sc.pp.log1p(g)
    panel = _load_panel()
    expression = g.X
    all_genes = list(g.var_names)
    gene_index = {gene: i for i, gene in enumerate(all_genes)}
    missing = sorted({gn for genes in panel.values() for gn in genes} - set(gene_index))
    print(f"panel genes missing from gene activity: {len(missing)}")
    usable = {ct: [gn for gn in genes if gn in gene_index] for ct, genes in panel.items()}
    usable = {ct: genes for ct, genes in usable.items() if genes}
    if not usable:
        analysis_notes.append("panel scoring skipped: no panel gene found in gene activity (check identifier type)")
    else:
        annotation_panel = list(panel.keys())
        gl = g.obs[group].astype(str).to_numpy() if group in g.obs else obs_col(adata, group).astype(str)
        annotation_scores = []
        for cluster in sorted_levels(gl):
            mask = gl == cluster
            units = int(mask.sum())
            if units == 0:
                continue
            mean = np.asarray(expression[mask].mean(axis=0)).ravel()
            total_mean = float(mean.sum())
            scores = []
            for cell_type, genes in panel.items():
                found = [gene_index[gn] for gn in genes if gn in gene_index]
                if not found:
                    scores.append({"cell_type": cell_type, "score": None, "genes_found": 0, "genes_missing": len(genes)})
                    continue
                panel_mean = float(mean[found].mean())
                rest_count = len(all_genes) - len(found)
                rest_mean = (total_mean - float(mean[found].sum())) / rest_count if rest_count else 0.0
                scores.append({"cell_type": cell_type, "score": round(panel_mean - rest_mean, 4),
                               "panel_mean": round(panel_mean, 4), "genes_found": len(found),
                               "genes_missing": len(genes) - len(found)})
            scored = sorted((s for s in scores if s["score"] is not None), key=lambda s: -s["score"])
            best = scored[0] if scored and scored[0]["score"] >= min_score else None
            annotation_scores.append({"cluster": cluster, "cells": units, "best": best, "scores": scores})
            labels_by_cluster[cluster] = best["cell_type"] if best else "unassigned"
        annotation_automated = "panel-score"

# 3) reference label transfer in gene-activity space
if not labels_by_cluster and reference_path and gene_adata is not None:
    if not os.path.exists(reference_path):
        analysis_notes.append("reference transfer skipped: reference not found: " + reference_path)
    else:
        try:
            from sklearn.neighbors import KNeighborsClassifier
            ref = sc.read_h5ad(reference_path)
            if reference_label_key not in ref.obs:
                raise ValueError(f"reference has no obs column '{reference_label_key}'")
            g = _gene_memory()
            shared = sorted(set(g.var_names) & set(ref.var_names))
            if len(shared) < 50:
                raise ValueError(f"only {len(shared)} shared genes; too few for transfer")
            gs = g[:, shared].copy()
            rs = ref[:, shared].copy()
            for obj in (gs, rs):
                sc.pp.normalize_total(obj, target_sum=1e4)
                sc.pp.log1p(obj)
            combined = sc.concat([rs, gs], join="inner")
            sc.pp.pca(combined, n_comps=min(30, len(shared) - 1))
            emb = combined.obsm["X_pca"]
            nref = rs.n_obs
            train_y = ref.obs[reference_label_key].astype(str).to_numpy()
            knn = KNeighborsClassifier(n_neighbors=min(15, nref), metric="cosine").fit(emb[:nref], train_y)
            pred = knn.predict(emb[nref:])
            gl = obs_col(adata, group).astype(str)
            for cluster in sorted_levels(gl):
                votes = pred[gl == cluster]
                labels_by_cluster[cluster] = str(pd.Series(votes).mode().iloc[0]) if len(votes) else "unassigned"
            annotation_automated = "reference-transfer"
            print(f"reference transfer on {len(shared)} shared genes from {nref} reference cells")
        except Exception as error:
            analysis_notes.append("reference transfer skipped: " + str(error))
            info("reference transfer skipped: " + str(error))

if not labels_by_cluster:
    analysis_notes.append(
        "annotation skipped: pass gene_anno (gene activity) for panel scoring, a reference_path, or a mapping_path"
    )
    info("annotation skipped: no gene activity, reference or mapping")
else:
    _cl = obs_col(adata, group).astype(str)
    _cell_labels = np.array([labels_by_cluster.get(c, "unassigned") for c in _cl])
    adata.obs[label_column] = _cell_labels
    annotation_assignments = [
        {"cluster": c, "label": l, "cells": int((_cl == c).sum())} for c, l in labels_by_cluster.items()
    ]
    analysis_steps.append("annotate")
    print(f"annotation ({annotation_automated}) -> obs['{label_column}']:")
    for entry in annotation_assignments:
        print(f"  cluster {entry['cluster']}: {entry['label']} ({entry['cells']} cells)")