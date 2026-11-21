# Batch integration removes technical variation between samples so clusters reflect
# biology, not which run a cell came from. The sample is the experimental unit:
# batch_column names each cell's sample. Harmony (default) corrects the spectral
# embedding; scanorama and mnc are alternatives. Integration needs at least two
# batches and a computed X_spectral, otherwise there is nothing to correct and the
# cell reports the gap. The corrected embedding (X_spectral_harmony / _scanorama /
# _mnn) is added to obsm; re-cluster on it (group_by / use_rep) if batch drove clusters.
if not obs_has(adata, batch_column):
    analysis_notes.append(f"integration skipped: no '{batch_column}' obs column (need a sample/batch label)")
    info(f"integration skipped: no {batch_column} column")
elif "X_spectral" not in embeddings:
    analysis_notes.append("integration skipped: no X_spectral; run the embed step first")
    info("integration skipped: no X_spectral")
else:
    _batches = obs_col(adata, batch_column).astype(str)
    _uniq = np.unique(_batches)
    if len(_uniq) < 2:
        analysis_notes.append(f"integration skipped: only one batch ({_uniq[0]}); nothing to integrate")
        info("integration skipped: single batch")
    else:
        try:
            if integrate_method == "harmony":
                snap.pp.harmony(adata, batch=batch_column, use_rep="X_spectral", inplace=True)
                _key = "X_spectral_harmony"
            elif integrate_method == "scanorama":
                snap.pp.scanorama_integrate(adata, batch=batch_column, use_rep="X_spectral", inplace=True)
                _key = "X_spectral_scanorama"
            elif integrate_method == "mnc":
                snap.pp.mnc_correct(adata, batch=batch_column, use_rep="X_spectral", inplace=True)
                _key = "X_spectral_mnn"
            else:
                raise ValueError("integrate_method must be harmony, scanorama or mnc")
            integration_info = {
                "method": integrate_method, "batch_column": batch_column,
                "n_batches": int(len(_uniq)), "batches": _uniq.tolist(), "embedding": _key,
            }
            embeddings.append(_key)
            analysis_steps.append("integrate")
            print(f"integrated {len(_uniq)} batches with {integrate_method} -> obsm['{_key}']")
        except Exception as error:
            analysis_notes.append("integration skipped: " + str(error))
            info("integration skipped: " + str(error))