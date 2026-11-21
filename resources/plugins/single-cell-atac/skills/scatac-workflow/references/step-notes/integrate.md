Purpose: align multiple samples so shared cell states cluster together despite
batch effects. Principle: Harmony (or Scanorama / MNN) corrects the spectral
embedding across the batch_column, producing a corrected embedding
(X_spectral_harmony etc.) for joint clustering and visualization. Interpretation:
compare cluster composition across samples before and after; integration should
mix batches within a cell type without erasing genuine condition differences.
Pitfalls: integration needs a batch_column with at least two samples, otherwise
the step reports a gap; corrected values are for clustering/visualization ONLY —
never feed them into peak calling or differential accessibility, which must return
to raw counts; over-correction can remove real biology, so do not integrate away
the condition you want to test; the sample is the experimental unit — fewer than
two biological replicates per condition means no condition-level claim. The
method (harmony/scanorama/mnc) is a choice; say which you used.