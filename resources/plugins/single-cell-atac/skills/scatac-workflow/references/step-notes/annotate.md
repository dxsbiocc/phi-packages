Purpose: turn each cluster from an ID into a cell-type hypothesis. Principle:
three routes, in priority order — an explicit cluster->label mapping; panel
scoring (panel mean minus the mean of everything else on the gene activity
matrix, so panel sizes stay comparable); or kNN label transfer from a scRNA
reference on shared genes. Interpretation: read the per-cluster scores and the
dot plot before accepting a call; a label needs positive markers AND absence of
competing-type markers; an unassigned cluster (below min_score) is information.
Pitfalls: a script must not invent labels — without mapping_path, panel_path or
reference_path nothing is annotated and the gap is reported; the built-in panel
is a generic example, real tissue needs a curated panel; gene activity is a proxy
so the label is inference not observation; label transfer needs a scRNA reference
(from the single-cell plugin) with enough shared genes, and a spatial spot mixes
cells but a scATAC barcode is one cell — still, accessibility is sparse so the
call is probabilistic. transfer_labels in snapatac2 is within-object kNN
propagation, not cross-modal ATAC<->RNA transfer; the cross-modal route here is
the explicit kNN classifier on shared-gene PCs.