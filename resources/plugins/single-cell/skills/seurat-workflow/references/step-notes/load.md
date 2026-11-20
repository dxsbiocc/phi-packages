Purpose: read the input and make its state explicit before any analysis.
Principle: a 10x directory or .rds becomes a Seurat object with raw counts in
the counts layer; the printed genes-per-cell distribution gates a min_genes
value that would delete every cell (Seurat's own error for that is opaque).
Interpretation: read assays/layers/shape — QC and pseudobulk need raw counts; an
object that arrives already normalized must be re-imported, not normalized
twice. Pitfalls: a filtered/subset matrix, Ensembl IDs instead of symbols, or a
shallow library all look like "bad quality" here; check the identifier type
before blaming the cells.
