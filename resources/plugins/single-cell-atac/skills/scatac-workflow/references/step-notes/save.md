Purpose: persist the analyzed object and release the backed files cleanly.
Principle: adata.write saves the object (the tile/peak matrix, obs cluster labels,
embeddings and uns results) to the output .h5ad; the separate peak_adata and
gene_adata objects are closed. Interpretation: the saved object is the handoff to
the next tool or to the single-cell plugin; the printed cells/features/steps/genome
are the audit summary. Pitfalls: a backed snapatac2 object MUST be closed before
the process exits or the HDF5 file corrupts — this step closes every backed
object it owns; the output path is inside the project, never the read-only plugin
environment; gene activity is a proxy and motif hits are predictions, so the
saved labels/annotations carry that caveat, not ground truth. Re-opening the
object later detects its kind and reports which steps are already done.