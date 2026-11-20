Purpose: estimate the cell-type mixture inside each spot/bin from a scRNA-seq
reference. Principle: a Visium spot (or a large Stereo-seq bin) captures several
cells, so its expression is a weighted sum of reference cell-type profiles;
deconvolution solves for those weights. The plugin uses a non-negative
least-squares fit against pseudobulk reference profiles (no deep-learning
dependency). Interpretation: read the per-unit composition table and map each
type's proportion in space; proportions are ESTIMATES, not observed cells.
Pitfalls: this REQUIRES a scRNA reference from the single-cell plugin — without it
there is no deconvolution; the reference's cell-type granularity and tissue match
bound the answer (a missing type cannot be recovered); small bins/cellbin are
near single-cell and deconvolution adds little; cell2location/tangram are
project-env gaps, not available here.