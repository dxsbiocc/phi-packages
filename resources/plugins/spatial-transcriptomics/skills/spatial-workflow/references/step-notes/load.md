Purpose: read the spatial input and make its platform, analysis unit and state
explicit before any analysis. Principle: a Visium Space Ranger directory/.h5, a
Stereo-seq GEF/GEM, a Visium HD bin table or an .h5ad becomes AnnData with raw
counts in the counts layer and tissue coordinates in obsm['spatial']; for
Stereo-seq the square bin aggregates N×N DNBs (bin_x = (x//N)*N) and cellbin
reads the cells SAW segmented upstream. Interpretation: read the platform/bin
report and the genes-per-unit distribution — it gates a min_genes value that
would delete every unit. Choosing a bin is a biological decision (large bins for
tissue domains, cellbin/small bins for cellular questions), not a default.
Pitfalls: an object without obsm['spatial'] supports only non-spatial steps; a
filtered matrix, Ensembl IDs instead of symbols, or a too-small bin all look like
"bad quality" here. Segmentation is never done in this cell — cellbin consumes
SAW cellCut output.