Purpose: reduce the sparse tile/peak matrix to a low-dimensional embedding for
clustering. Principle: snapatac2.tl.spectral computes a diffusion-map / spectral
embedding on the feature graph — the scATAC analogue of LSI/PCA in scRNA, suited
to binary-ish accessibility data; a kNN graph is then built on X_spectral.
Interpretation: n_comps controls the embedding dimensionality (default 30);
read the spectral eigenvalues to see how many components carry structure.
Pitfalls: this is spectral embedding, not TF-IDF+LSI literally — do not claim
"LSI" was run; the embedding is on selected features, so feature selection
upstream changes it; corrected/integrated embeddings are for visualization only
and never feed peak calling or differential tests. The kNN graph (n_neighbors)
defines the neighbourhood Leiden will cluster — a different k changes the
clusters.