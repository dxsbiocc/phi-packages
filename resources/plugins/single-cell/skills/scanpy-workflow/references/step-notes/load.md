Purpose: read the input and make its state explicit before any analysis.
Principle: a 10x directory, .h5, .h5ad, .loom or count table becomes AnnData with
raw counts snapshotted into the counts layer; the printed genes-per-cell
distribution gates a min_genes value that would delete every cell.
Interpretation: read the shape and layers — QC, HVG (seurat_v3) and pseudobulk
need raw counts; an object that arrives already normalized must be re-imported,
not normalized twice. Pitfalls: a filtered/subset matrix, Ensembl IDs instead of
symbols, or a shallow library all look like "bad quality" here; check the
identifier type before blaming the cells.
