Purpose: flag two-cells-in-one-droplet artifacts before they invent clusters.
Principle: scDblFinder simulates artificial doublets from the observed counts
and scores each cell against them; heterotypic doublets (two different cell
types) masquerade as novel "hybrid" states, homotypic ones are nearly invisible.
Interpretation: the expected rate scales with loaded cells (~4% at ~9k loaded);
inspect the scores on the UMAP and run detection per sample when libraries were
captured separately. Pitfalls: scDblFinder is not in the plugin environment, so
this cell degrades to an explicit message — an unremoved doublet is information
you need before trusting the clusters; the rate is a loading estimate, not a
free knob.
