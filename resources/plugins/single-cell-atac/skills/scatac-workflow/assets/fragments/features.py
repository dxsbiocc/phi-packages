# The feature matrix for clustering. Fragments import with no feature matrix, so
# bin the genome into uniform tiles (bin_size bp) — tiles cover the genome evenly
# and are robust for the initial embedding. A 10x peak/tile matrix already carries
# features, so skip tiling there. select_features then drops uninformative bins
# (too rare or too ubiquitous) and removes ENCODE blacklist artifact regions, which
# otherwise dominate the embedding with technical signal.
if adata.n_vars == 0:
    snap.pp.add_tile_matrix(adata, bin_size=bin_size, inplace=True, n_jobs=n_jobs)
    analysis_steps.append("tile-matrix")
    print(f"tile matrix: {adata.n_vars} tiles at {bin_size} bp")
else:
    info(f"input already has {adata.n_vars} features; skipping tile matrix")

_blacklist = resolve_blacklist()
if _blacklist:
    info("blacklist: " + _blacklist)
else:
    analysis_notes.append("no blacklist applied (genome has no inline ENCODE list and no blacklist path given)")
snap.pp.select_features(
    adata, n_features=min(n_features, adata.n_vars), blacklist=_blacklist,
    inplace=True, n_jobs=n_jobs,
)
analysis_steps.append("select-features")
_sel = var_col(adata, "selected")
_n_sel = int(np.asarray(_sel).sum()) if _sel is not None else adata.n_vars
print(f"selected {_n_sel} features for the embedding")

# Doublet removal: scrublet simulates doublets from the selected features and scores
# each cell; filter_doublets drops cells above the probability threshold. This needs
# the feature matrix, which is why it lives here and not in the threshold-only filter
# step. expected_doublet_rate scales with how many cells were loaded.
if doublet_filter:
    _before_dbl = adata.n_obs
    try:
        snap.pp.scrublet(
            adata, features="selected", expected_doublet_rate=expected_doublet_rate,
            random_state=seed, inplace=True,
        )
        snap.pp.filter_doublets(adata, probability_threshold=doublet_threshold, inplace=True)
        _removed = _before_dbl - adata.n_obs
        qc_filters.append({"filter": f"scrublet doublet > {doublet_threshold}", "removed": int(_removed)})
        analysis_steps.append("doublet-filter")
        print(f"scrublet removed {_removed} predicted doublets of {_before_dbl} cells")
        if adata.n_obs == 0:
            raise ValueError("scrublet removed every cell; raise doublet_threshold or disable doublet_filter")
    except ValueError:
        raise
    except Exception as error:
        analysis_notes.append("scrublet skipped: " + str(error))
        info("scrublet skipped: " + str(error))
else:
    info("doublet_filter disabled; skipping scrublet")