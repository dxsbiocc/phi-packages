Purpose: compress the variable-gene space into orthogonal axes ordered by
variance. Principle: PCA on the scaled matrix with the highly-variable mask
separates biological signal (first axes) from technical noise (later axes) while
the object keeps all genes so markers remain testable. Interpretation: read the
elbow and the first-ten variance share, then keep a few PCs beyond the bend so
rare populations survive; dims in the config cell is that decision. Pitfalls:
too few PCs merges rare types into their neighbours; too many re-admits noise
into the neighbour graph.
