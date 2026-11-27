# TF regulon activity inference (AUCell recovery AUC and weighted ULM linear model
# t-statistic on the unscaled lognorm layer). Expression-based regulon activity is a
# transcription factor activity proxy from target co-expression, not physical DNA
# binding (pair with scATAC / ChIP-seq to confirm motif occupancy).
import matplotlib.pyplot as plt
from scipy import stats as scipy_stats


def _bh_adjust(pvals):
    arr = np.asarray(pvals, dtype=np.float64)
    n = len(arr)
    if n == 0:
        return arr
    order = np.argsort(arr)
    ranked = arr[order] * n / np.arange(1, n + 1, dtype=np.float64)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty_like(ranked)
    out[order] = np.clip(ranked, 0.0, 1.0)
    return out


def _compute_aucell_scores(expr_dense, target_indices_list, top_frac=0.05):
    """Compute normalized AUCell area-under-recovery-curve scores per cell."""
    n_cells, n_genes = expr_dense.shape
    n_top = max(5, min(n_genes, int(np.ceil(float(top_frac) * n_genes))))
    scores = np.zeros((n_cells, len(target_indices_list)), dtype=np.float64)
    batch = 512
    for start in range(0, n_cells, batch):
        end = min(n_cells, start + batch)
        block = expr_dense[start:end]
        # Descending rank per cell (0 = highest expressed gene)
        desc_order = np.argsort(-block, axis=1, kind="mergesort")
        ranks = np.empty_like(desc_order)
        rows = np.arange(end - start)[:, None]
        ranks[rows, desc_order] = np.arange(n_genes)[None, :]
        for j, idx_arr in enumerate(target_indices_list):
            k = len(idx_arr)
            max_auc = float(min(k, n_top) * (2 * n_top - min(k, n_top) + 1) / 2.0)
            if max_auc <= 0:
                continue
            r_sub = ranks[:, idx_arr]
            # Contribution of each target ranked within [0, n_top): (n_top - r) summed
            contrib = np.maximum(0.0, float(n_top) - r_sub)
            scores[start:end, j] = contrib.sum(axis=1) / max_auc
    return scores


def _compute_ulm_scores(expr_dense, weight_matrix):
    """Vectorized weighted univariate linear model t-statistics per cell x regulon."""
    n_cells, n_genes = expr_dense.shape
    if n_genes < 3:
        return np.zeros((n_cells, weight_matrix.shape[1]), dtype=np.float64)
    w_mean = weight_matrix.mean(axis=0, keepdims=True)
    w_c = weight_matrix - w_mean
    denom_w = np.sum(w_c ** 2, axis=0)
    y_mean = expr_dense.mean(axis=1, keepdims=True)
    y_c = expr_dense - y_mean
    denom_y = np.sum(y_c ** 2, axis=1, keepdims=True)
    cov = y_c @ w_c
    df = max(1, n_genes - 2)
    scores = np.zeros((n_cells, weight_matrix.shape[1]), dtype=np.float64)
    for j in range(weight_matrix.shape[1]):
        dw = denom_w[j]
        if dw <= 1e-12:
            continue
        r2_num = (cov[:, j:j + 1] ** 2)
        rss = np.maximum(1e-12, denom_y * dw - r2_num)
        t_val = cov[:, j:j + 1] * np.sqrt(df / rss)
        scores[:, j] = t_val[:, 0]
    return scores


if regulon_method not in ("aucell", "ulm"):
    raise ValueError(f"regulon_method must be 'aucell' or 'ulm', got: {regulon_method}")

_curated_regulons = {
    "TP53": ["CDKN1A", "MDM2", "BAX", "BBC3", "GADD45A", "FAS", "TNFRSF10B", "DDB2", "RRM2B", "TIGAR"],
    "NFKB1": ["NFKBIA", "TNFAIP3", "IL1B", "IL6", "CCL2", "CXCL8", "ICAM1", "BCL2A1", "TRAF1", "RELB"],
    "STAT1": ["IRF1", "GBP1", "CXCL9", "CXCL10", "ISG15", "STAT2", "TAP1", "PSMB9", "IDO1", "OAS1"],
    "STAT3": ["SOCS3", "BCL3", "MCL1", "MYC", "CCND1", "VEGFA", "BIRC5", "PIM1", "IL10", "BCL2L1"],
    "MYC": ["ODC1", "NCL", "NPM1", "CDK4", "LDHA", "CAD", "SRM", "FBL", "EIF4E", "APEX1"],
    "HIF1A": ["VEGFA", "SLC2A1", "LDHA", "PDK1", "BNIP3", "CA9", "ENO1", "HK2", "PGK1", "ALDOA"],
    "JUN": ["FOS", "ATF3", "DUSP1", "EGR1", "CDKN1A", "MMP1", "IL6", "CCL2", "HB-EGF", "CCND1"],
    "FOS": ["JUN", "EGR1", "NR4A1", "DUSP1", "ATF3", "PTGS2", "IL6", "MMP1", "VEGFD", "IER2"],
    "IRF7": ["ISG15", "IFIT1", "IFIT3", "MX1", "OAS1", "STAT1", "IFI44L", "RSAD2", "IRF1", "HERC5"],
    "E2F1": ["MCM2", "MCM6", "PCNA", "TOP2A", "MKI67", "CCNE1", "CDK1", "RRM2", "TYMS", "CDC6"],
    "FOXP3": ["IL2RA", "CTLA4", "TIGIT", "TNFRSF18", "IKZF2", "CCR8", "LRRC32", "ICOS", "ENTPD1", "SELL"],
    "TBX21": ["IFNG", "GZMB", "PRF1", "CXCR3", "CCL5", "EOMES", "NKG7", "GNLY", "KLRD1", "STAT4"],
    "GATA3": ["IL4", "IL5", "IL13", "CCR4", "MAF", "IL1RL1", "PTGDR2", "KLRG1", "AREG", "BATF"],
    "PAX5": ["CD19", "MS4A1", "CD79A", "CD79B", "BLNK", "EBF1", "SPIB", "BANK1", "FCRLA", "VPREB3"],
    "SPI1": ["CSF1R", "TYROBP", "FCER1G", "CD68", "CD14", "LYZ", "AIF1", "ITGAM", "CLEC7A", "CTSS"],
    "CEBPB": ["IL6", "CCL3", "PTGS2", "S100A8", "S100A9", "SOD2", "NAMPT", "CD14", "FCGR3A", "CXCL2"],
    "SOX9": ["COL2A1", "ACAN", "COMP", "FOXA2", "KRT19", "SOX4", "SFRP5", "FGFBP1", "SPP1", "CLU"],
    "KLF4": ["CDKN1A", "CDH1", "EPCAM", "KRT8", "KRT18", "IVL", "SPRR1A", "CEACAM1", "MUC1", "ZBTB7B"],
}

_reg_dict = {}
if regulons_path:
    if not os.path.exists(regulons_path):
        raise ValueError("regulons_path not found: " + regulons_path)
    if regulons_path.lower().endswith(".gmt"):
        with open(regulons_path, encoding="utf-8") as _fh:
            for _line in _fh:
                _parts = [p.strip() for p in _line.rstrip("\r\n").split("\t") if p.strip()]
                if len(_parts) >= 3:
                    _reg_dict[_parts[0]] = [(g, 1.0) for g in _parts[2:]]
    else:
        _sep = "," if regulons_path.lower().endswith(".csv") else "\t"
        _rdf = pd.read_csv(regulons_path, sep=_sep)
        _rcols = {str(c).strip().lower(): c for c in _rdf.columns}
        _tf_col = _rcols.get("source") or _rcols.get("tf") or _rcols.get("regulon")
        _tgt_col = _rcols.get("target") or _rcols.get("gene")
        _wt_col = _rcols.get("weight") or _rcols.get("mor")
        if not _tf_col or not _tgt_col:
            raise ValueError("regulons_path table requires source/tf and target/gene columns")
        for _, _row in _rdf.iterrows():
            _tf = str(_row[_tf_col]).strip()
            _tgt = str(_row[_tgt_col]).strip()
            _w = float(_row[_wt_col]) if (_wt_col and pd.notna(_row[_wt_col])) else 1.0
            if _tf and _tgt:
                _reg_dict.setdefault(_tf, []).append((_tgt, _w))
    _reg_resource = os.path.basename(regulons_path)
else:
    _reg_dict = {tf: [(g, 1.0) for g in genes] for tf, genes in _curated_regulons.items()}
    _reg_resource = "curated_tf_starter_v1"

_var_lookup = {}
for _idx, _vn in enumerate(adata.var_names):
    _var_lookup.setdefault(str(_vn).upper(), _idx)

_active_tfs = []
_target_indices_list = []
_weight_rows = []
for _tf, _pairs in sorted(_reg_dict.items()):
    _matched = {}
    for _g, _w in _pairs:
        _gi = _var_lookup.get(str(_g).upper())
        if _gi is not None:
            _matched[_gi] = float(_w)
    if len(_matched) >= int(regulon_min_targets):
        _active_tfs.append(_tf)
        _target_indices_list.append(sorted(_matched.keys()))
        _w_col = np.zeros(adata.n_vars, dtype=np.float64)
        for _gi, _w in _matched.items():
            _w_col[_gi] = _w
        _weight_rows.append(_w_col)

if not _active_tfs:
    info("regulon skipped: 0 TFs met regulon_min_targets=" + str(regulon_min_targets))
    analysis_notes.append(
        f"regulon skipped: no TFs had >= {regulon_min_targets} targets matched in var_names"
    )
    reg_info = {
        "skipped": "no regulons matched min_targets",
        "method": regulon_method,
        "resource": _reg_resource,
    }
else:
    _mat_src = adata.layers["lognorm"] if "lognorm" in adata.layers else adata.X
    _dense = np.asarray(_mat_src.toarray() if hasattr(_mat_src, "toarray") else _mat_src, dtype=np.float64)
    if regulon_method == "aucell":
        _scores = _compute_aucell_scores(_dense, _target_indices_list, top_frac=regulon_top_frac)
    else:
        _w_mat = np.column_stack(_weight_rows)
        _scores = _compute_ulm_scores(_dense, _w_mat)

    adata.obsm["X_regulon"] = _scores
    adata.uns["regulon_names"] = list(_active_tfs)
    _score_df = pd.DataFrame(_scores, index=adata.obs_names, columns=_active_tfs)

    _reg_dir = os.path.dirname(os.path.abspath(regulon_output))
    if _reg_dir:
        os.makedirs(_reg_dir, exist_ok=True)
    _score_df.to_csv(regulon_output, sep="\t")

    _group_col = group_by if group_by in adata.obs.columns else None
    _summary_rows = []
    if _group_col:
        _gvals = adata.obs[_group_col].astype(str)
        _glevels = sorted_levels(_gvals)
        _tf_global_mean = _scores.mean(axis=0)
        _tf_global_std = np.where(_scores.std(axis=0) > 1e-12, _scores.std(axis=0), 1.0)
        _raw_pvals = []
        for _lvl in _glevels:
            _in_mask = (_gvals == _lvl).to_numpy()
            _out_mask = ~_in_mask
            for _ti, _tf in enumerate(_active_tfs):
                _v_in = _scores[_in_mask, _ti]
                _v_out = _scores[_out_mask, _ti]
                _mean_in = float(_v_in.mean())
                _mean_out = float(_v_out.mean()) if _out_mask.any() else _mean_in
                if len(_v_in) >= 2 and len(_v_out) >= 2:
                    _, _pval = scipy_stats.mannwhitneyu(_v_in, _v_out, alternative="two-sided")
                else:
                    _pval = 1.0
                _raw_pvals.append(float(_pval))
                _summary_rows.append({
                    "group": _lvl,
                    "tf": _tf,
                    "n_cells": int(_in_mask.sum()),
                    "n_targets": len(_target_indices_list[_ti]),
                    "mean_activity": round(_mean_in, 5),
                    "diff_vs_rest": round(_mean_in - _mean_out, 5),
                    "z_score": round((_mean_in - float(_tf_global_mean[_ti])) / float(_tf_global_std[_ti]), 4),
                    "pvalue": float(_pval),
                })
        _fdrs = _bh_adjust(_raw_pvals)
        for _r, _q in zip(_summary_rows, _fdrs):
            _r["fdr"] = float(_q)

        _sum_df = pd.DataFrame(_summary_rows).sort_values(["fdr", "diff_vs_rest"], ascending=[True, False])
        _summary_path = os.path.join(_reg_dir or ".", "regulon_group_summary.tsv")
        _sum_df.to_csv(_summary_path, sep="\t", index=False)

        _piv = _sum_df.pivot(index="tf", columns="group", values="z_score").reindex(_active_tfs).fillna(0.0)
        _top_tfs = _piv.var(axis=1).sort_values(ascending=False).head(min(25, len(_active_tfs))).index.tolist()
        _plot_df = _piv.loc[_top_tfs]
        _fig, _ax = plt.subplots(figsize=(max(5.5, 0.7 * _plot_df.shape[1] + 2.0), max(4.2, 0.35 * len(_top_tfs) + 1.5)))
        _im = _ax.imshow(_plot_df.to_numpy(), aspect="auto", cmap="RdBu_r", vmin=-2.0, vmax=2.0)
        _ax.set_xticks(np.arange(_plot_df.shape[1]))
        _ax.set_xticklabels(_plot_df.columns.tolist(), rotation=30, ha="right", fontsize=8)
        _ax.set_yticks(np.arange(len(_top_tfs)))
        _ax.set_yticklabels(_top_tfs, fontsize=8)
        _ax.set_title(f"TF regulon z-score ({regulon_method.upper()})")
        _fig.colorbar(_im, ax=_ax, label="Group mean z-score")
        show_plot(_ax, "regulon-heatmap")
        _top_records = _sum_df.head(15).to_dict(orient="records")
    else:
        _summary_path = None
        _top_records = []

    analysis_steps.append("regulon")
    reg_info = {
        "method": regulon_method,
        "resource": _reg_resource,
        "n_regulons": len(_active_tfs),
        "n_cells": int(adata.n_obs),
        "group_by": _group_col,
        "regulons_tested": _active_tfs,
        "top_regulons": _top_records,
        "tables": [p for p in (regulon_output, _summary_path) if p],
    }
    print(f"regulon ({regulon_method}, {len(_active_tfs)} TFs across {adata.n_obs} cells)")
