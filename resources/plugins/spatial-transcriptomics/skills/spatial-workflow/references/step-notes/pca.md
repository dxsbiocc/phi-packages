Purpose: reduce the variable-gene space to orthogonal axes ordered by variance.
Principle: PCA rotates the HVG space so the first components capture the most
variation; the mask keeps all genes in the object (markers and SVG need them)
while the axes are computed from variable genes only. Interpretation: read the
elbow — axes past the bend are mostly noise, but keep a few extra so rare domains
survive; the first-ten-percent summary tells you how much structure the leading
axes carry. Pitfalls: too few PCs merges distinct domains, too many adds noise
that Leiden then clusters on; dims (PCs used downstream) should not exceed the
computed n_pcs. PCA is in expression space, never on tissue coordinates.