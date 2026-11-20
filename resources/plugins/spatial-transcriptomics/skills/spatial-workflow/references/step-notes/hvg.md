Purpose: keep the genes that carry signal and drop the noise before PCA.
Principle: most genes vary only with technical noise; seurat_v3 models the
mean-variance trend of raw counts directly, so it runs on the counts layer (not
the log matrix) and ranks genes by excess variability. Interpretation: these
HVGs drive PCA/clustering; the lognorm snapshot keeps ALL genes so markers and
SVG can still test any gene. Spatially variable genes (the svg step) are a
separate, coordinate-aware selection — HVGs vary across units, SVGs vary across
space. Pitfalls: too few HVGs erases rare domains, too many dilutes the axes;
seurat_v3 needs scikit-misc (shipped here). Read the mean-variance plot before
trusting n_hvg.