# Cis-regulatory Peak-to-Gene linkage correlating peak accessibility against
# paired scRNA expression within peak2gene_max_dist bp of gene TSSs, with
# trans-chromosome empirical null calibration and BH FDR correction.
import gzip
import os
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


def _parse_region(region_str):
    s = str(region_str).strip()
    if ":" in s and "-" in s:
        chrom, span = s.split(":", 1)
        parts = span.split("-", 1)
        try:
            a, b = int(parts[0]), int(parts[1])
            return chrom, a, b, (a + b) // 2
        except ValueError:
            return None
    parts = s.rsplit("-", 2)
    if len(parts) == 3:
        try:
            a, b = int(parts[1]), int(parts[2])
            return parts[0], a, b, (a + b) // 2
        except ValueError:
            return None
    return None


def _load_tss_table(tss_path=None, gtf_path=None):
    rows = {}
    if tss_path and os.path.exists(tss_path):
        sep = "," if tss_path.lower().endswith(".csv") else "\t"
        df = pd.read_csv(tss_path, sep=sep)
        cols = {str(c).strip().lower(): c for c in df.columns}
        g_col = cols.get("gene") or cols.get("gene_name") or cols.get("symbol") or df.columns[0]
        c_col = cols.get("chrom") or cols.get("chr") or df.columns[1]
        t_col = cols.get("tss") or cols.get("start") or df.columns[2]
        for _, r in df.iterrows():
            rows[str(r[g_col]).strip()] = (str(r[c_col]).strip(), int(r[t_col]))
        return rows
    if gtf_path and os.path.exists(gtf_path):
        opener = gzip.open if gtf_path.lower().endswith(".gz") else open
        with opener(gtf_path, "rt", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\r\n").split("\t")
                if len(parts) < 9 or parts[2] != "gene":
                    continue
                chrom, start_s, end_s, strand, attr = parts[0], parts[3], parts[4], parts[6], parts[8]
                tss = int(end_s) if strand == "-" else int(start_s)
                gname = None
                for seg in attr.split(";"):
                    seg = seg.strip()
                    if seg.startswith("gene_name"):
                        gname = seg.split(None, 1)[-1].strip('"\' ')
                        break
                    if seg.startswith("gene_id") and gname is None:
                        gname = seg.split(None, 1)[-1].strip('"\' ')
                if gname:
                    rows.setdefault(gname, (chrom, tss))
    return rows


def _read_matrix_any(path):
    if path.lower().endswith(".h5ad"):
        import anndata as ad
        obj = ad.read_h5ad(path)
        x = np.asarray(obj.X.toarray() if hasattr(obj.X, "toarray") else obj.X, dtype=np.float64)
        return [str(v) for v in obj.var_names], [str(v) for v in obj.obs_names], x, obj.obs
    sep = "," if path.lower().endswith((".csv", ".csv.gz")) else "\t"
    df = pd.read_csv(path, sep=sep, index_col=0)
    return [str(x) for x in df.index], [str(x) for x in df.columns], df.to_numpy(dtype=np.float64).T, None


def _compute_peak2gene_links(atac_mat, peak_names, atac_obs_names, rna_mat, gene_names, rna_obs_names,
                              tss_map, max_dist=250000, n_null=200, seed=42,
                              atac_groups=None, rna_groups=None):
    """Compute cis peak-gene Pearson/Spearman correlations and trans empirical FDR."""
    shared_cells = [c for c in atac_obs_names if c in set(rna_obs_names)]
    if len(shared_cells) >= 5:
        a_idx = {c: i for i, c in enumerate(atac_obs_names)}
        r_idx = {c: i for i, c in enumerate(rna_obs_names)}
        A = atac_mat[[a_idx[c] for c in shared_cells], :]
        R = rna_mat[[r_idx[c] for c in shared_cells], :]
        pairing_mode = "paired_barcodes"
        n_units = len(shared_cells)
    elif atac_groups is not None and rna_groups is not None:
        shared_grps = sorted(set(atac_groups) & set(rna_groups))
        if len(shared_grps) < 3:
            raise ValueError("peak2gene requires >= 5 shared cell barcodes or >= 3 shared groups")
        A = np.vstack([atac_mat[atac_groups == g].mean(axis=0) for g in shared_grps])
        R = np.vstack([rna_mat[rna_groups == g].mean(axis=0) for g in shared_grps])
        pairing_mode = "group_pseudobulk"
        n_units = len(shared_grps)
    else:
        raise ValueError("No shared cell barcodes or shared group labels between ATAC and RNA matrices")

    parsed_peaks = {}
    for pi, pk in enumerate(peak_names):
        coords = _parse_region(pk)
        if coords is not None:
            parsed_peaks[pi] = coords

    gene_idx_by_name = {g: gi for gi, g in enumerate(gene_names)}
    cand_pairs = []
    for gname, (gchrom, gtss) in tss_map.items():
        gi = gene_idx_by_name.get(gname)
        if gi is None or R[:, gi].std() <= 1e-12:
            continue
        for pi, (pchrom, pstart, pend, pmid) in parsed_peaks.items():
            if pchrom != gchrom or A[:, pi].std() <= 1e-12:
                continue
            dist = int(pmid - gtss)
            if abs(dist) <= int(max_dist):
                cand_pairs.append((pi, gi, gname, peak_names[pi], pchrom, dist))

    if not cand_pairs:
        return pd.DataFrame(), pairing_mode, n_units

    rng = np.random.default_rng(int(seed))
    peak_indices = sorted(parsed_peaks.keys())
    null_corrs = []
    genes_active = list({c[1] for c in cand_pairs})
    while len(null_corrs) < int(n_null) and peak_indices and genes_active:
        pi_r = int(rng.choice(peak_indices))
        gi_r = int(rng.choice(genes_active))
        if parsed_peaks[pi_r][0] == tss_map.get(gene_names[gi_r], ("", 0))[0]:
            if len({parsed_peaks[k][0] for k in peak_indices}) > 1:
                continue
        if A[:, pi_r].std() <= 1e-12 or R[:, gi_r].std() <= 1e-12:
            continue
        r_null, _ = scipy_stats.pearsonr(A[:, pi_r], R[:, gi_r])
        if np.isfinite(r_null):
            null_corrs.append(float(r_null))
    null_arr = np.asarray(null_corrs, dtype=np.float64)

    rows = []
    for pi, gi, gname, pname, chrom, dist in cand_pairs:
        r_val, p_par = scipy_stats.pearsonr(A[:, pi], R[:, gi])
        rho_val, _ = scipy_stats.spearmanr(A[:, pi], R[:, gi])
        if len(null_arr) >= 10:
            p_emp = float((np.sum(np.abs(null_arr) >= abs(r_val)) + 1.0) / (len(null_arr) + 1.0))
        else:
            p_emp = float(p_par)
        rows.append({
            "peak": pname,
            "gene": gname,
            "chrom": chrom,
            "tss_distance_bp": int(dist),
            "abs_distance_bp": int(abs(dist)),
            "is_promoter": abs(dist) <= 2000,
            "pearson_r": round(float(r_val), 5),
            "spearman_rho": round(float(rho_val), 5),
            "pvalue": p_emp,
        })

    out_df = pd.DataFrame(rows)
    out_df["fdr"] = _bh_adjust(out_df["pvalue"].to_numpy())
    out_df = out_df.sort_values(["fdr", "pearson_r"], ascending=[True, False])
    return out_df, pairing_mode, n_units


_src_adata = peak_adata if peak_adata is not None else adata
if _src_adata is not None:
    _peak_names = [str(v) for v in _src_adata.var_names]
    _atac_cells = [str(v) for v in _src_adata.obs_names]
    _ax = _src_adata.X
    _atac_mat = np.asarray(_ax.toarray() if hasattr(_ax, "toarray") else _ax, dtype=np.float64)
    _atac_grps = obs_col(_src_adata, group_by)
elif data_path and os.path.exists(data_path):
    _peak_names, _atac_cells, _atac_mat, _ = (
        _read_matrix_any(data_path)
    )
    _atac_grps = None
    if meta_path and os.path.exists(meta_path):
        _msep = "," if meta_path.lower().endswith(".csv") else "\t"
        _mdf = pd.read_csv(meta_path, sep=_msep, index_col=0)
        if group_by in _mdf.columns:
            _atac_grps = _mdf[group_by].reindex(_atac_cells).fillna("unknown").astype(str).to_numpy()
else:
    raise ValueError("peak2gene requires an ATAC peak AnnData object or peak x cell TSV/CSV")

if not rna_input or not os.path.exists(rna_input):
    analysis_notes.append("peak2gene skipped: rna_input matrix is missing")
    peak2gene_info = {"skipped": "missing rna_input"}
else:
    _gene_names, _rna_cells, _rna_mat, _rna_obs = _read_matrix_any(rna_input)
    _rna_grps = None
    if _rna_obs is not None and group_by in _rna_obs.columns:
        _rna_grps = _rna_obs[group_by].astype(str).to_numpy()
    elif rna_meta and os.path.exists(rna_meta):
        _rm_sep = "," if rna_meta.lower().endswith(".csv") else "\t"
        _rm_df = pd.read_csv(rna_meta, sep=_rm_sep, index_col=0)
        if group_by in _rm_df.columns:
            _rna_grps = _rm_df[group_by].reindex(_rna_cells).fillna("unknown").astype(str).to_numpy()

    _tss_map = _load_tss_table(tss_path=gene_tss, gtf_path=gene_anno)
    if not _tss_map:
        analysis_notes.append("peak2gene skipped: provide gene_tss TSV/BED or gene_anno GTF/GFF3")
        peak2gene_info = {"skipped": "missing gene_tss or gene_anno"}
    else:
        _links_df, _pairing_mode, _n_units = _compute_peak2gene_links(
            _atac_mat, _peak_names, _atac_cells,
            _rna_mat, _gene_names, _rna_cells,
            _tss_map,
            max_dist=peak2gene_max_dist,
            n_null=peak2gene_n_null,
            seed=seed,
            atac_groups=_atac_grps,
            rna_groups=_rna_grps,
        )
        _out_dir = os.path.dirname(os.path.abspath(peak2gene_output))
        if _out_dir:
            os.makedirs(_out_dir, exist_ok=True)
        _links_df.to_csv(peak2gene_output, sep="\t", index=False)

        _n_sig = int((_links_df["fdr"] <= 0.05).sum()) if len(_links_df) else 0
        if figures_dir and len(_links_df) > 0:
            os.makedirs(figures_dir, exist_ok=True)
            _fig, _axes = plt.subplots(1, 2, figsize=(10.5, 4.2))
            _sig_mask = (_links_df["fdr"] <= 0.05).to_numpy()
            _axes[0].scatter(
                _links_df["tss_distance_bp"].to_numpy() / 1000.0,
                _links_df["pearson_r"].to_numpy(),
                c=np.where(_sig_mask, "#E11D48", "#94A3B8"),
                s=16, alpha=0.75,
            )
            _axes[0].axhline(0, color="#475569", lw=0.8, ls="--")
            _axes[0].set_xlabel("Distance to TSS (kb)")
            _axes[0].set_ylabel("Peak-Gene Pearson r")
            _axes[0].set_title(f"Cis Peak2Gene links (sig={_n_sig}/{len(_links_df)})")

            _top = _links_df.head(min(15, len(_links_df)))
            _lbls = [f"{r['gene']} | {r['peak']}" for _, r in _top.iterrows()]
            _axes[1].barh(np.arange(len(_top)), _top["pearson_r"].to_numpy()[::-1], color="#0EA5E9")
            _axes[1].set_yticks(np.arange(len(_top)))
            _axes[1].set_yticklabels(_lbls[::-1], fontsize=7)
            _axes[1].set_xlabel("Pearson r")
            _axes[1].set_title("Top cis links")
            _fig.tight_layout()
            _fig_path = os.path.join(figures_dir, "peak2gene-summary.png")
            _fig.savefig(_fig_path, dpi=150, bbox_inches="tight")
            plt.close(_fig)
            figure_paths["peak2gene-summary"] = _fig_path

        analysis_steps.append("peak2gene")
        peak2gene_info = {
            "pairing_mode": _pairing_mode,
            "n_paired_units": int(_n_units),
            "max_distance_bp": int(peak2gene_max_dist),
            "n_links_tested": int(len(_links_df)),
            "n_significant_links": _n_sig,
            "n_promoter_links": int(_links_df["is_promoter"].sum()) if len(_links_df) else 0,
            "output_links": peak2gene_output,
            "top_links": _links_df.head(15).to_dict(orient="records") if len(_links_df) else [],
        }
        print(f"peak2gene ({_pairing_mode}, n={_n_units}): {_n_sig}/{len(_links_df)} significant links")
