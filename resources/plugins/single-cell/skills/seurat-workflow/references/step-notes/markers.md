Purpose: produce the evidence a cell-type call must rest on. Principle:
FindAllMarkers tests each cluster against the rest (Wilcoxon by default, the
best-calibrated rank test); pct.1/pct.2 say how specific a gene is, not just how
different. Interpretation: read top genes with their pct columns and the
heatmap; downsample unequal clusters (max.cells.per.ident) so a huge cluster
does not dominate the test. Pitfalls: clusters were defined from the same data,
so p-values are double-dipped and descriptive only — rank, do not claim
significance; tests run on log-normalized full genes, never on integrated or
scaled values; avg_log2FC is this comparison's fold change, not a portable
effect size.
