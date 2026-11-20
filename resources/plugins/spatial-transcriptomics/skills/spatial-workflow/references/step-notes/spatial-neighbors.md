Purpose: build the tissue-coordinate adjacency graph every spatial statistic
reads. Principle: two units can be neighbours in tissue even when their
expression differs, so this graph is built on obsm[coord_key], NOT on the
expression PCA/UMAP. neighbor_method picks the geometry: grid (Visium hex/square
lattice, n_neighs=6), knn (k nearest units; Stereo-seq bins or irregular
layouts), radius (everything within a physical distance), delaunay
(triangulation, adaptive). Interpretation: read the edge count and mean degree —
a graph that is too sparse misses real neighbours, too dense bridges gaps.
Pitfalls: the adjacency definition changes downstream conclusions, so state the
method and parameters; grid assumes a regular lattice (prefer knn/radius for
irregular bins); a null radius is derived as 1.5× the median nearest-neighbour
spacing. squidpy 1.8 deprecates the catch-all spatial_neighbors, so the
mode-specific builder is called directly.