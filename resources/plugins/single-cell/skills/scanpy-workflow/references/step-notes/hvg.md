Purpose: keep the genes that carry structure and drop the noise that dilutes
every distance. Principle: seurat_v3 models raw counts directly (hence the
counts layer) and ranks genes by variance beyond the mean-variance trend of
technical noise. Interpretation: read the top genes and the mean-variance plot —
ribosomal, mitochondrial or one dominant population monopolizing the list is a
warning; with batches, per-batch ranking stops one sample choosing every gene.
Pitfalls: HVG on scaled or log values with a counts flavor is a category error;
too few HVGs merges rare types, too many re-admits noise.
