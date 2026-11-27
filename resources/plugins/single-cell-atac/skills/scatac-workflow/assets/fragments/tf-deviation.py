# GC- and accessibility-matched chromVAR-style per-cell TF motif deviations
# (Schep et al. 2017): expected accessibility E = n_i * e_j, 5x5 quantile binning
# on (GC, log10(mean_access)), K matched background peak draws, raw deviation
# (O - E) / E, and background mean/SD bias-corrected z-scores.
import gzip
import os
import matplotlib.pyplot as plt
from scipy import stats as scipy_stats


def _open_maybe_gz(path):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace") if path.lower().endswith(".gz") else open(path, "r", encoding="utf-8", errors="replace")


def _parse_peak_table(path):
    """Read a peak x cell TSV/CSV table -> (peaks, cells, X_cell_by_peak float64, obs_df_or_None)."""
    sep = "," if path.lower().endswith((".csv", ".csv.gz")) else "\t"
    df = pd.read_csv(path, sep=sep, index_col=0)
    peaks = [str(x) for x in df.index.tolist()]
    cells = [str(x) for x in df.columns.tolist()]
    return peaks, cells, df.to_numpy(dtype=np.float64).T


def _read_fasta_peaks(fasta_path, peak_names):
    seqs = {}
    chrom_seqs = {}
    cur_name = None
    buf = []
    with _open_maybe_gz(fasta_path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur_name is not None:
                    chrom_seqs[cur_name] = "".join(buf).upper()
                cur_name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
        if cur_name is not None:
            chrom_seqs[cur_name] = "".join(buf).upper()

    for pk in peak_names:
        if pk in chrom_seqs:
            seqs[pk] = chrom_seqs[pk]
            continue
        if ":" in pk and "-" in pk:
            chrom, span = pk.split(":", 1)
            s_str, e_str = span.split("-", 1)
            try:
                s, e = max(0, int(s_str)), int(e_str)
            except ValueError:
                continue
            cseq = chrom_seqs.get(chrom) or chrom_seqs.get(chrom.lstrip("chr"))
            if cseq is not None and e <= len(cseq):
                seqs[pk] = cseq[s:e]
    return seqs


def _parse_meme_pwms(meme_path):
    motifs = {}
    cur_name = None
    rows = []
    in_mat = False
    with open(meme_path, "r", encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            line = raw.strip()
            if line.startswith("MOTIF "):
                if cur_name and rows:
                    motifs[cur_name] = np.asarray(rows, dtype=np.float64)
                parts = line.split(None, 2)
                cur_name = parts[2] if len(parts) > 2 else parts[1]
                rows = []
                in_mat = False
            elif line.startswith("letter-probability matrix:"):
                in_mat = True
                rows = []
            elif in_mat:
                if not line or line.startswith(("URL", "MOTIF")):
                    in_mat = False
                    continue
                vals = line.split()
                if len(vals) == 4:
                    try:
                        rows.append([float(v) for v in vals])
                    except ValueError:
                        in_mat = False
                else:
                    in_mat = False
    if cur_name and rows:
        motifs[cur_name] = np.asarray(rows, dtype=np.float64)
    return motifs


def _scan_pwm_matches(peak_names, peak_seqs, motifs, score_threshold=6.0):
    base_idx = {"A": 0, "C": 1, "G": 2, "T": 3}
    comp = {0: 3, 1: 2, 2: 1, 3: 0}
    tf_names = sorted(motifs.keys())
    matches = np.zeros((len(peak_names), len(tf_names)), dtype=np.float64)
    log_odds_list = []
    for tf in tf_names:
        pwm = np.clip(motifs[tf], 1e-4, 1.0)
        pwm = pwm / pwm.sum(axis=1, keepdims=True)
        lo_f = np.log2(pwm / 0.25)
        lo_r = lo_f[::-1, [comp[i] for i in range(4)]]
        log_odds_list.append((lo_f, lo_r))

    for pi, pk in enumerate(peak_names):
        seq = peak_seqs.get(pk, "")
        if not seq:
            continue
        arr = np.array([base_idx.get(ch, -1) for ch in seq], dtype=np.int32)
        for ti, (lo_f, lo_r) in enumerate(log_odds_list):
            w = lo_f.shape[0]
            if len(arr) < w:
                continue
            best = -1e9
            for start in range(len(arr) - w + 1):
                win = arr[start:start + w]
                if (win < 0).any():
                    continue
                sf = float(lo_f[np.arange(w), win].sum())
                sr = float(lo_r[np.arange(w), win].sum())
                if max(sf, sr) > best:
                    best = max(sf, sr)
            if best >= score_threshold:
                matches[pi, ti] = 1.0
    return tf_names, matches


def _chromvar_deviations(counts_cell_x_peak, peak_motif_mat, peak_gc=None, n_background=50, seed=42):
    """Exact GC- and accessibility-matched chromVAR z-scores and variability."""
    X = np.asarray(counts_cell_x_peak, dtype=np.float64)
    M = (np.asarray(peak_motif_mat, dtype=np.float64) > 0).astype(np.float64)
    n_cells, n_peaks = X.shape
    if n_cells < 2 or n_peaks < 5:
        raise ValueError("chromVAR requires at least 2 cells and 5 peaks")

    total_fragments_per_cell = X.sum(axis=1, keepdims=True)
    total_fragments_per_cell = np.where(total_fragments_per_cell > 0, total_fragments_per_cell, 1.0)
    sum_total = float(total_fragments_per_cell.sum())
    peak_sum = X.sum(axis=0, keepdims=True)
    peak_frac = peak_sum / max(1.0, sum_total)  # 1 x n_peaks
    expected_per_cell = total_fragments_per_cell @ peak_frac  # n_cells x n_peaks

    mean_access = X.mean(axis=0)
    log_acc = np.log10(mean_access + 1e-6)
    acc_bins = pd.qcut(pd.Series(log_acc).rank(method="first"), q=min(5, max(1, n_peaks // 5)), labels=False).to_numpy()
    acc_only_pools = {}
    for idx, ab in enumerate(acc_bins):
        acc_only_pools.setdefault(int(ab), []).append(idx)

    if peak_gc is not None:
        gc_bins = pd.qcut(pd.Series(peak_gc).rank(method="first"), q=min(5, max(1, n_peaks // 5)), labels=False).to_numpy()
        bin_ids = acc_bins * 10 + gc_bins
    else:
        bin_ids = acc_bins

    raw_pools = {}
    for idx, bid in enumerate(bin_ids):
        raw_pools.setdefault(int(bid), []).append(idx)
    bin_to_indices = {}
    for idx, bid in enumerate(bin_ids):
        p = raw_pools[int(bid)]
        if len(p) >= 3:
            bin_to_indices[int(bid)] = p
        else:
            acc_p = acc_only_pools.get(int(acc_bins[idx]), [])
            bin_to_indices[int(bid)] = acc_p if len(acc_p) >= 3 else list(range(n_peaks))
    all_indices = np.arange(n_peaks)

    rng = np.random.default_rng(int(seed))
    K = max(10, int(n_background))
    n_tfs = M.shape[1]
    z_mat = np.zeros((n_cells, n_tfs), dtype=np.float64)
    raw_dev_mat = np.zeros((n_cells, n_tfs), dtype=np.float64)

    for ti in range(n_tfs):
        fg_idx = np.where(M[:, ti] > 0)[0]
        if len(fg_idx) == 0:
            continue
        obs_m = X[:, fg_idx].sum(axis=1)
        exp_m = np.maximum(1e-12, expected_per_cell[:, fg_idx].sum(axis=1))
        raw_d = (obs_m - exp_m) / exp_m
        raw_dev_mat[:, ti] = raw_d

        bg_devs = np.zeros((n_cells, K), dtype=np.float64)
        for k_iter in range(K):
            bg_idx = np.empty(len(fg_idx), dtype=np.int64)
            for p_i, pk in enumerate(fg_idx):
                pool = bin_to_indices.get(int(bin_ids[pk]), all_indices)
                bg_idx[p_i] = int(rng.choice(pool))
            obs_b = X[:, bg_idx].sum(axis=1)
            exp_b = np.maximum(1e-12, expected_per_cell[:, bg_idx].sum(axis=1))
            bg_devs[:, k_iter] = (obs_b - exp_b) / exp_b

        bg_mean = bg_devs.mean(axis=1)
        bg_sd = bg_devs.std(axis=1, ddof=1)
        bias_corr = raw_d - bg_mean
        z_mat[:, ti] = np.where(bg_sd > 1e-12, bias_corr / bg_sd, 0.0)

    variability = z_mat.std(axis=0, ddof=1) if n_cells > 1 else np.zeros(n_tfs)
    return z_mat, raw_dev_mat, variability


_src_adata = peak_adata if peak_adata is not None else adata
_peak_names = None
_cell_names = None
_counts_mat = None
_cell_groups = None

if _src_adata is not None:
    _peak_names = [str(v) for v in _src_adata.var_names]
    _cell_names = [str(v) for v in _src_adata.obs_names]
    _raw_x = _src_adata.X
    _counts_mat = np.asarray(_raw_x.toarray() if hasattr(_raw_x, "toarray") else _raw_x, dtype=np.float64)
    try:
        _gvals = obs_col(_src_adata, group_by)
        if _gvals is not None:
            _cell_groups = np.asarray([str(x) for x in _gvals])
    except Exception:
        _cell_groups = None
elif data_path and os.path.exists(data_path):
    _peak_names, _cell_names, _counts_mat = _parse_peak_table(data_path)
    if meta_path and os.path.exists(meta_path):
        _msep = "," if meta_path.lower().endswith(".csv") else "\t"
        _mdf = pd.read_csv(meta_path, sep=_msep, index_col=0)
        if group_by in _mdf.columns:
            _cell_groups = _mdf[group_by].reindex(_cell_names).fillna("unknown").astype(str).to_numpy()
else:
    raise ValueError("tf-deviation requires a peak AnnData object or a peak x cell TSV/CSV matrix")

_peak_gc_arr = None
_gc_mode = "accessibility_only_fallback"
if peak_gc and os.path.exists(peak_gc):
    _gsep = "," if peak_gc.lower().endswith(".csv") else "\t"
    _gc_df = pd.read_csv(peak_gc, sep=_gsep)
    _gc_cols = {str(c).strip().lower(): c for c in _gc_df.columns}
    _pk_c = _gc_cols.get("peak") or _gc_df.columns[0]
    _gc_c = _gc_cols.get("gc") or _gc_df.columns[1]
    _gc_map = dict(zip(_gc_df[_pk_c].astype(str), _gc_df[_gc_c].astype(float)))
    _peak_gc_arr = np.array([float(_gc_map.get(pk, 0.5)) for pk in _peak_names], dtype=np.float64)
    _gc_mode = "peak_gc_table"
elif genome_fasta and os.path.exists(genome_fasta):
    _peak_seqs_gc = _read_fasta_peaks(genome_fasta, _peak_names)
    if _peak_seqs_gc:
        _peak_gc_arr = np.array([
            ((s.count("G") + s.count("C")) / max(1, len(s))) if (s := _peak_seqs_gc.get(pk)) else 0.5
            for pk in _peak_names
        ], dtype=np.float64)
        _gc_mode = "genome_fasta_gc"

if _gc_mode == "accessibility_only_fallback":
    analysis_notes.append(
        "tf-deviation: no peak_gc table or genome_fasta provided; background matched on accessibility bins only"
    )

_tf_names = []
_motif_mat = None
if motif_matches and os.path.exists(motif_matches):
    _msep = "," if motif_matches.lower().endswith((".csv", ".csv.gz")) else "\t"
    _mm_df = pd.read_csv(motif_matches, sep=_msep, index_col=0)
    _mm_df = _mm_df.reindex(_peak_names).fillna(0.0)
    _tf_names = [str(c) for c in _mm_df.columns]
    _motif_mat = _mm_df.to_numpy(dtype=np.float64)
elif genome_fasta and os.path.exists(genome_fasta) and motif_db and os.path.exists(motif_db):
    _seqs = _read_fasta_peaks(genome_fasta, _peak_names)
    _pwms = _parse_meme_pwms(motif_db)
    if not _seqs or not _pwms:
        analysis_notes.append("tf-deviation skipped: empty peak sequences or MEME motifs")
        tfdev_info = {"skipped": "empty sequences or motifs"}
    else:
        _tf_names, _motif_mat = _scan_pwm_matches(_peak_names, _seqs, _pwms)
else:
    analysis_notes.append("tf-deviation skipped: provide motif_matches TSV/CSV or genome_fasta + motif_db (.meme)")
    tfdev_info = {"skipped": "missing motif_matches or genome_fasta + motif_db"}

if _motif_mat is not None and len(_tf_names) > 0:
    _keep_tfs = [ti for ti in range(len(_tf_names)) if (_motif_mat[:, ti] > 0).sum() >= 1]
    if not _keep_tfs:
        analysis_notes.append("tf-deviation skipped: 0 TFs matched any peaks")
        tfdev_info = {"skipped": "no TF matched peaks"}
    else:
        _tf_names = [_tf_names[i] for i in _keep_tfs]
        _motif_mat = _motif_mat[:, _keep_tfs]
        _z_mat, _raw_mat, _var_vec = _chromvar_deviations(
            _counts_mat, _motif_mat, peak_gc=_peak_gc_arr, n_background=n_background, seed=seed
        )
        if _src_adata is not None:
            _src_adata.obsm["X_tf_deviation"] = _z_mat
            _src_adata.uns["tf_deviation_names"] = list(_tf_names)

        _out_dir = os.path.dirname(os.path.abspath(tfdev_output))
        if _out_dir:
            os.makedirs(_out_dir, exist_ok=True)
        _dev_df = pd.DataFrame(_z_mat, index=_cell_names, columns=_tf_names)
        _dev_df.to_csv(tfdev_output, sep="\t")

        _var_df = pd.DataFrame({
            "tf": _tf_names,
            "n_peaks_matched": [int((_motif_mat[:, i] > 0).sum()) for i in range(len(_tf_names))],
            "variability_sd": np.round(_var_vec, 5),
            "mean_abs_z": np.round(np.abs(_z_mat).mean(axis=0), 5),
            "max_z": np.round(_z_mat.max(axis=0), 5),
        }).sort_values("variability_sd", ascending=False)
        _var_path = os.path.join(_out_dir or ".", "tf_variability.tsv")
        _var_df.to_csv(_var_path, sep="\t", index=False)

        _fig_path = None
        if figures_dir:
            os.makedirs(figures_dir, exist_ok=True)
            _fig, _axes = plt.subplots(1, 2, figsize=(10.5, 4.2))
            _top_v = _var_df.head(min(15, len(_var_df)))
            _axes[0].barh(np.arange(len(_top_v)), _top_v["variability_sd"].to_numpy()[::-1], color="#0284c7")
            _axes[0].set_yticks(np.arange(len(_top_v)))
            _axes[0].set_yticklabels(_top_v["tf"].tolist()[::-1], fontsize=8)
            _axes[0].set_xlabel("chromVAR deviation SD across cells")
            _axes[0].set_title("Top variable TF motifs")

            _top_tfs = _var_df["tf"].head(min(15, len(_var_df))).tolist()
            if _cell_groups is not None and len(set(_cell_groups)) >= 2:
                _grp_means = _dev_df[_top_tfs].groupby(_cell_groups).mean().T
                _im = _axes[1].imshow(_grp_means.to_numpy(), aspect="auto", cmap="RdBu_r", vmin=-3, vmax=3)
                _axes[1].set_xticks(np.arange(_grp_means.shape[1]))
                _axes[1].set_xticklabels(_grp_means.columns.tolist(), rotation=30, ha="right", fontsize=8)
                _axes[1].set_yticks(np.arange(len(_top_tfs)))
                _axes[1].set_yticklabels(_top_tfs, fontsize=8)
                _axes[1].set_title(f"Mean motif deviation z ({group_by})")
                _fig.colorbar(_im, ax=_axes[1], label="Mean z-score")
            else:
                _axes[1].hist(_var_df["variability_sd"].to_numpy(), bins=min(20, max(5, len(_var_df))), color="#0f766e")
                _axes[1].set_xlabel("Deviation SD")
                _axes[1].set_title("TF variability distribution")
            _fig.tight_layout()
            _fig_path = os.path.join(figures_dir, "tf-deviation-summary.png")
            _fig.savefig(_fig_path, dpi=150, bbox_inches="tight")
            plt.close(_fig)
            figure_paths["tf-deviation-summary"] = _fig_path

        analysis_steps.append("tf-deviation")
        tfdev_info = {
            "method": "chromvar_gc_matched",
            "gc_matching": _gc_mode,
            "n_background": int(n_background),
            "n_cells": int(_z_mat.shape[0]),
            "n_peaks": int(_counts_mat.shape[1]),
            "n_tfs": int(len(_tf_names)),
            "output_deviations": tfdev_output,
            "output_variability": _var_path,
            "top_variable_tfs": _var_df.head(10).to_dict(orient="records"),
        }
        print(f"tf-deviation ({len(_tf_names)} TFs, {n_background} bg draws, {_gc_mode})")
