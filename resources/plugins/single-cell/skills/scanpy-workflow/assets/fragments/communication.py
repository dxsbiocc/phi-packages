# Cell-cell communication (CellPhoneDB-style heteromeric subunit-minimum permutation
# test on the unscaled lognorm layer). Expression co-occurrence suggests potential
# signaling capacity between cell populations, NOT physical contact (which requires
# spatial transcriptomics or imaging validation).
import matplotlib.pyplot as plt

if group_by not in adata.obs.columns:
    raise ValueError(
        f"group column '{group_by}' not in obs; available: " + ", ".join(adata.obs.columns)
    )

_groups_series = adata.obs[group_by].astype(str)
_levels = sorted_levels(_groups_series)
if len(_levels) < 2:
    info("communication needs >= 2 cell groups; skipping")
    analysis_notes.append("communication skipped: fewer than 2 groups in " + group_by)
    ccc_info = {"skipped": "fewer than 2 groups", "group_by": group_by}
else:
    if interactions_path:
        if not os.path.exists(interactions_path):
            raise ValueError("interactions_path not found: " + interactions_path)
        _sep = "," if interactions_path.lower().endswith(".csv") else "\t"
        _lr_df = pd.read_csv(interactions_path, sep=_sep)
        _cols_low = {str(c).strip().lower(): c for c in _lr_df.columns}
        if "ligand" in _cols_low and "receptor" in _cols_low:
            _raw_pairs = list(
                zip(_lr_df[_cols_low["ligand"]].astype(str), _lr_df[_cols_low["receptor"]].astype(str))
            )
        elif "source" in _cols_low and "target" in _cols_low:
            _raw_pairs = list(
                zip(_lr_df[_cols_low["source"]].astype(str), _lr_df[_cols_low["target"]].astype(str))
            )
        else:
            raise ValueError("interactions_path requires ligand/receptor or source/target columns")
        _resource_name = os.path.basename(interactions_path)
    else:
        # Curated starter resource covering major cytokines, chemokines, growth factors,
        # ECM-integrin heterodimers, Notch/Wnt, and immune checkpoints (subunit '+' syntax).
        _raw_pairs = [
            ("TGFB1", "TGFBR1+TGFBR2"), ("TGFB2", "TGFBR1+TGFBR2"), ("BMP4", "BMPR1A+BMPR2"),
            ("VEGFA", "FLT1"), ("VEGFA", "KDR"), ("PDGFB", "PDGFRB"), ("EGF", "EGFR"),
            ("AREG", "EGFR"), ("HGF", "MET"), ("IGF1", "IGF1R"), ("FGF2", "FGFR1"),
            ("GAS6", "AXL"), ("MIF", "CD74+CXCR4"), ("MIF", "CD74+CD44"), ("APP", "CD74"),
            ("TNF", "TNFRSF1A"), ("TNF", "TNFRSF1B"), ("IL1B", "IL1R1+IL1RAP"),
            ("IL6", "IL6R+IL6ST"), ("IL10", "IL10RA+IL10RB"), ("IFNG", "IFNGR1+IFNGR2"),
            ("IL2", "IL2RA+IL2RB+IL2RG"), ("IL15", "IL15RA+IL2RB+IL2RG"), ("IL4", "IL4R"),
            ("IL7", "IL7R+IL2RG"), ("IL33", "IL1RL1"), ("CSF1", "CSF1R"), ("FLT3LG", "FLT3"),
            ("CCL2", "CCR2"), ("CCL3", "CCR1"), ("CCL4", "CCR5"), ("CCL5", "CCR5"),
            ("CXCL8", "CXCR1"), ("CXCL8", "CXCR2"), ("CXCL9", "CXCR3"), ("CXCL10", "CXCR3"),
            ("CXCL12", "CXCR4"), ("CXCL13", "CXCR5"), ("CCL19", "CCR7"), ("CCL21", "CCR7"),
            ("CX3CL1", "CX3CR1"), ("XCL1", "XCR1"),
            ("FN1", "ITGAV+ITGB1"), ("FN1", "ITGA5+ITGB1"), ("COL1A1", "ITGA1+ITGB1"),
            ("COL1A2", "ITGA2+ITGB1"), ("SPP1", "CD44"), ("SPP1", "ITGAV+ITGB1"),
            ("ICAM1", "ITGAL+ITGB2"), ("VCAM1", "ITGA4+ITGB1"), ("THBS1", "CD47"),
            ("LGALS9", "HAVCR2"), ("CD274", "PDCD1"), ("PDCD1LG2", "PDCD1"),
            ("CD80", "CD28"), ("CD86", "CD28"), ("CD80", "CTLA4"), ("CD86", "CTLA4"),
            ("CD40LG", "CD40"), ("FASLG", "FAS"), ("TNFSF10", "TNFRSF10B"),
            ("HLA-E", "KLRC1+KLRD1"), ("NECTIN2", "TIGIT"), ("PVR", "TIGIT"),
            ("DLL1", "NOTCH1"), ("JAG1", "NOTCH1"), ("JAG1", "NOTCH2"),
            ("WNT5A", "FZD5"), ("ANXA1", "FPR1"), ("C3", "C3AR1"), ("C5", "C5AR1"),
        ]
        _resource_name = "curated_starter_v1"

    # Case-insensitive mapping to adata.var_names so human and mouse symbols both match
    _var_lookup = {}
    for _vn in adata.var_names:
        _var_lookup.setdefault(str(_vn).upper(), str(_vn))

    _resolved_pairs = []
    for _lig_raw, _rec_raw in _raw_pairs:
        _l_subs = [s.strip() for s in str(_lig_raw).replace("_", "+").split("+") if s.strip()]
        _r_subs = [s.strip() for s in str(_rec_raw).replace("_", "+").split("+") if s.strip()]
        if not _l_subs or not _r_subs:
            continue
        if all(s.upper() in _var_lookup for s in _l_subs) and all(s.upper() in _var_lookup for s in _r_subs):
            _l_matched = [_var_lookup[s.upper()] for s in _l_subs]
            _r_matched = [_var_lookup[s.upper()] for s in _r_subs]
            _resolved_pairs.append(("+".join(_l_matched), "+".join(_r_matched), _l_matched, _r_matched))

    if not _resolved_pairs:
        info("communication skipped: 0 ligand-receptor pairs matched adata.var_names")
        analysis_notes.append(
            "communication skipped: no ligand-receptor subunits matched var_names (check gene symbols vs Ensembl IDs or supply interactions_path)"
        )
        ccc_info = {
            "skipped": "no ligand-receptor pairs matched var_names",
            "resource": _resource_name,
            "group_by": group_by,
        }
    else:
        _unique_genes = sorted({g for _, _, ls, rs in _resolved_pairs for g in (ls + rs)})
        _gene_to_col = {g: i for i, g in enumerate(_unique_genes)}
        _src_layer = adata.layers["lognorm"] if "lognorm" in adata.layers else adata.X
        _sub_view = adata[:, _unique_genes]
        _sub_mat = _sub_view.layers["lognorm"] if "lognorm" in _sub_view.layers else _sub_view.X
        if hasattr(_sub_mat, "toarray"):
            _X_lr = np.asarray(_sub_mat.toarray(), dtype=np.float64)
        else:
            _X_lr = np.asarray(_sub_mat, dtype=np.float64)
        _X_pos = (_X_lr > 0).astype(np.float64)

        _n_cells = adata.n_obs
        _n_grps = len(_levels)
        _grp_to_idx = {g: i for i, g in enumerate(_levels)}
        _cell_codes = np.array([_grp_to_idx[g] for g in _groups_series], dtype=np.int64)
        _grp_counts = np.bincount(_cell_codes, minlength=_n_grps).astype(np.float64)
        _safe_counts = np.maximum(_grp_counts, 1.0)[:, None]

        # Observed group means and expressed fractions (G x M)
        _one_hot = np.zeros((_n_cells, _n_grps), dtype=np.float64)
        _one_hot[np.arange(_n_cells), _cell_codes] = 1.0
        _obs_means = (_one_hot.T @ _X_lr) / _safe_counts
        _obs_pcts = (_one_hot.T @ _X_pos) / _safe_counts

        _n_pairs = len(_resolved_pairs)
        _lig_mean_obs = np.zeros((_n_grps, _n_pairs), dtype=np.float64)
        _rec_mean_obs = np.zeros((_n_grps, _n_pairs), dtype=np.float64)
        _lig_ok = np.zeros((_n_grps, _n_pairs), dtype=bool)
        _rec_ok = np.zeros((_n_grps, _n_pairs), dtype=bool)

        for _pi, (_, _, _ls, _rs) in enumerate(_resolved_pairs):
            _li = [_gene_to_col[g] for g in _ls]
            _ri = [_gene_to_col[g] for g in _rs]
            _lig_mean_obs[:, _pi] = np.min(_obs_means[:, _li], axis=1)
            _rec_mean_obs[:, _pi] = np.min(_obs_means[:, _ri], axis=1)
            _lig_ok[:, _pi] = np.min(_obs_pcts[:, _li], axis=1) >= ccc_min_pct
            _rec_ok[:, _pi] = np.min(_obs_pcts[:, _ri], axis=1) >= ccc_min_pct

        # Score tensor: sender A x receiver B x pair P
        _active_mask = _lig_ok[:, None, :] & _rec_ok[None, :, :]
        _score_obs = np.where(
            _active_mask,
            0.5 * (_lig_mean_obs[:, None, :] + _rec_mean_obs[None, :, :]),
            0.0,
        )

        # Vectorized cell-label permutation test
        _rng = np.random.default_rng(seed)
        _exceed = np.zeros_like(_score_obs, dtype=np.int32)
        _n_perms_eff = max(50, int(ccc_perms))
        for _ in range(_n_perms_eff):
            _perm_codes = _rng.permutation(_cell_codes)
            _oh_p = np.zeros((_n_cells, _n_grps), dtype=np.float64)
            _oh_p[np.arange(_n_cells), _perm_codes] = 1.0
            _p_means = (_oh_p.T @ _X_lr) / _safe_counts
            _lm_p = np.empty((_n_grps, _n_pairs), dtype=np.float64)
            _rm_p = np.empty((_n_grps, _n_pairs), dtype=np.float64)
            for _pi, (_, _, _ls, _rs) in enumerate(_resolved_pairs):
                _lm_p[:, _pi] = np.min(_p_means[:, [_gene_to_col[g] for g in _ls]], axis=1)
                _rm_p[:, _pi] = np.min(_p_means[:, [_gene_to_col[g] for g in _rs]], axis=1)
            _sc_p = 0.5 * (_lm_p[:, None, :] + _rm_p[None, :, :])
            _exceed += (_active_mask & (_sc_p >= _score_obs - 1e-12)).astype(np.int32)

        _pvals = np.where(_active_mask, (_exceed + 1.0) / (_n_perms_eff + 1.0), 1.0)

        # BH FDR over active tested interactions
        _flat_active = _active_mask.ravel()
        _fdrs_flat = np.ones(_pvals.size, dtype=np.float64)
        if np.any(_flat_active):
            _pv_act = _pvals.ravel()[_flat_active]
            _ord = np.argsort(_pv_act)
            _n_act = len(_pv_act)
            _adj = _pv_act[_ord] * _n_act / np.arange(1, _n_act + 1, dtype=np.float64)
            _adj = np.minimum.accumulate(_adj[::-1])[::-1]
            _adj = np.clip(_adj, 0.0, 1.0)
            _res_q = np.empty(_n_act, dtype=np.float64)
            _res_q[_ord] = _adj
            _fdrs_flat[_flat_active] = _res_q
        _fdrs = _fdrs_flat.reshape(_pvals.shape)

        _rows = []
        for _si, _sname in enumerate(_levels):
            for _ti, _tname in enumerate(_levels):
                for _pi, (_lname, _rname, _, _) in enumerate(_resolved_pairs):
                    if not _active_mask[_si, _ti, _pi]:
                        continue
                    _rows.append(
                        {
                            "source": _sname,
                            "target": _tname,
                            "ligand": _lname,
                            "receptor": _rname,
                            "interaction": f"{_lname} -> {_rname}",
                            "ligand_mean": round(float(_lig_mean_obs[_si, _pi]), 4),
                            "receptor_mean": round(float(_rec_mean_obs[_ti, _pi]), 4),
                            "score": round(float(_score_obs[_si, _ti, _pi]), 4),
                            "pvalue": round(float(_pvals[_si, _ti, _pi]), 5),
                            "fdr": round(float(_fdrs[_si, _ti, _pi]), 5),
                        }
                    )

        _ccc_df = pd.DataFrame(_rows)
        if not _ccc_df.empty:
            _ccc_df = _ccc_df.sort_values(["fdr", "pvalue", "score"], ascending=[True, True, False]).reset_index(drop=True)
        _sig_df = _ccc_df[_ccc_df["fdr"] <= ccc_fdr] if not _ccc_df.empty else _ccc_df

        _saved_table = None
        if ccc_output:
            os.makedirs(os.path.dirname(os.path.abspath(ccc_output)) or ".", exist_ok=True)
            _ccc_df.to_csv(ccc_output, sep="\t", index=False)
            _saved_table = ccc_output

        adata.uns["phi_communication"] = {
            "group_by": group_by,
            "resource": _resource_name,
            "min_pct": float(ccc_min_pct),
            "n_perms": int(_n_perms_eff),
            "fdr_threshold": float(ccc_fdr),
            "n_pairs_matched": int(_n_pairs),
            "n_active_tests": int(len(_ccc_df)),
            "n_significant": int(len(_sig_df)),
        }
        analysis_steps.append("communication")

        # Figure 1: Sender x Receiver significant interaction count heatmap
        _cnt_mat = np.zeros((_n_grps, _n_grps), dtype=int)
        if not _sig_df.empty:
            for _s, _t in zip(_sig_df["source"], _sig_df["target"]):
                _cnt_mat[_grp_to_idx[_s], _grp_to_idx[_t]] += 1

        _fig_hm, _ax_hm = plt.subplots(figsize=(max(4.6, 0.55 * _n_grps + 2.0), max(4.0, 0.48 * _n_grps + 1.6)))
        _im = _ax_hm.imshow(_cnt_mat, cmap="YlOrRd", aspect="auto")
        _ax_hm.set_xticks(np.arange(_n_grps))
        _ax_hm.set_xticklabels(_levels, rotation=40, ha="right", fontsize=8.5)
        _ax_hm.set_yticks(np.arange(_n_grps))
        _ax_hm.set_yticklabels(_levels, fontsize=8.5)
        for _i in range(_n_grps):
            for _j in range(_n_grps):
                _ax_hm.text(_j, _i, str(_cnt_mat[_i, _j]), ha="center", va="center", fontsize=8)
        _ax_hm.set_xlabel("Receiver population")
        _ax_hm.set_ylabel("Sender population")
        _ax_hm.set_title(f"Significant L–R Interactions (FDR ≤ {ccc_fdr})")
        _fig_hm.colorbar(_im, ax=_ax_hm, fraction=0.04, pad=0.03)
        show_plot(_ax_hm, "communication-heatmap")

        # Figure 2: Dotplot of top interactions
        _plot_df = (_sig_df if not _sig_df.empty else _ccc_df).head(25).copy()
        if not _plot_df.empty:
            _plot_df["pair_axis"] = _plot_df["source"] + " → " + _plot_df["target"]
            _y_labels = list(dict.fromkeys(_plot_df["interaction"].tolist()))[::-1]
            _x_labels = list(dict.fromkeys(_plot_df["pair_axis"].tolist()))
            _yi = {k: i for i, k in enumerate(_y_labels)}
            _xi = {k: i for i, k in enumerate(_x_labels)}

            _fig_dp, _ax_dp = plt.subplots(
                figsize=(max(5.5, 0.45 * len(_x_labels) + 2.6), max(4.2, 0.34 * len(_y_labels) + 1.8))
            )
            _xs = [_xi[v] for v in _plot_df["pair_axis"]]
            _ys = [_yi[v] for v in _plot_df["interaction"]]
            _sizes = 35.0 + 55.0 * (-np.log10(np.clip(_plot_df["fdr"].to_numpy(dtype=float), 1e-4, 1.0)))
            _sc = _ax_dp.scatter(
                _xs, _ys, s=_sizes, c=_plot_df["score"].to_numpy(dtype=float), cmap="viridis", alpha=0.9, edgecolors="#0F172A", linewidths=0.5
            )
            _ax_dp.set_xticks(np.arange(len(_x_labels)))
            _ax_dp.set_xticklabels(_x_labels, rotation=45, ha="right", fontsize=8)
            _ax_dp.set_yticks(np.arange(len(_y_labels)))
            _ax_dp.set_yticklabels(_y_labels, fontsize=8)
            _ax_dp.set_title("Top Ligand–Receptor Pairs (size = −log10 FDR, color = score)")
            _fig_dp.colorbar(_sc, ax=_ax_dp, fraction=0.035, pad=0.03, label="Mean L–R lognorm")
            show_plot(_ax_dp, "communication-dotplot")

        ccc_info = {
            "method": "heteromeric_permutation",
            "group_by": group_by,
            "resource": _resource_name,
            "min_pct": float(ccc_min_pct),
            "n_perms": int(_n_perms_eff),
            "fdr_threshold": float(ccc_fdr),
            "n_pairs_in_resource_matched": int(_n_pairs),
            "n_active_interactions_tested": int(len(_ccc_df)),
            "n_significant": int(len(_sig_df)),
            "interactions_table": _saved_table,
            "top_interactions": (_sig_df if not _sig_df.empty else _ccc_df).head(15).to_dict(orient="records"),
        }
        print(
            f"communication tested {len(_ccc_df)} active group-pair L-R combinations across {_n_pairs} matched pairs; "
            f"{len(_sig_df)} passed FDR <= {ccc_fdr}"
        )
