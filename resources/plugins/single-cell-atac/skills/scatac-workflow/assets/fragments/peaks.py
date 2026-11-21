# Per-cluster peak calling turns the uniform tiles into the open-chromatin regions
# that actually differ between cell states. MACS3 calls peaks on each cluster's
# fragments (it needs the fragments, so a 10x peak-matrix input without fragments
# reports a gap here). merge_peaks extends each summit by peak_half_width and merges
# overlaps into one non-overlapping set — the features of record for differential
# accessibility. The peak-count matrix is a SEPARATE object (peak_adata): the backed
# tile matrix cannot be resized in place, and keeping tiles for the embedding while
# peaks drive markers is the cleaner separation. obs (incl. the cluster labels) is
# copied onto the peak object so markers group by the same clusters.
group = group_by if obs_has(adata, group_by) else (cluster_columns[0] if cluster_columns else None)
if group is None:
    analysis_notes.append("peaks skipped: no cluster column; run the cluster step first")
    info("peaks skipped: no cluster column")
else:
    _blacklist = resolve_blacklist()
    try:
        snap.tl.macs3(
            adata, groupby=group, qvalue=peak_qvalue, call_broad_peaks=call_broad_peaks,
            blacklist=_blacklist, inplace=True, n_jobs=n_jobs,
        )
        _tables = adata.uns["macs3"] if "macs3" in adata.uns else {}
        peak_tables = {str(k): int(len(v)) for k, v in _tables.items()}
        sizes, _ = resolve_chrom_sizes()
        merged = snap.tl.merge_peaks(_tables, sizes, half_width=peak_half_width)
        merged_peaks = int(merged.height)
        if merged_peaks == 0:
            analysis_notes.append("no peaks called; relax peak_qvalue or check cells per cluster")
            info("no peaks called on this data")
        else:
            os.makedirs(tables_dir, exist_ok=True)
            bed = os.path.join(tables_dir, "merged_peaks.bed")
            with open(bed, "w", encoding="utf-8") as handle:
                for _p in merged["Peaks"].to_list():
                    _chrom, _rng = _p.split(":")
                    _s, _e = _rng.split("-")
                    handle.write(f"{_chrom}\t{_s}\t{_e}\n")
            _peak_path = os.path.join(tables_dir, "peak_matrix.h5ad")
            peak_adata = snap.pp.make_peak_matrix(
                adata, use_rep=merged["Peaks"].to_list(), inplace=False, file=_peak_path,
            )
            analysis_steps.append("peaks")
            print(f"peaks per cluster: {peak_tables}")
            print(f"merged {merged_peaks} non-overlapping peaks -> {_peak_path} (BED: {bed})")
    except Exception as error:
        analysis_notes.append("peak calling skipped: " + str(error))
        info("peak calling skipped: " + str(error))