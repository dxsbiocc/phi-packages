Purpose: map expression-space clusters and their markers back onto tissue — the
payoff of spatial data. Principle: the map uses obsm[coord_key] directly
(sc.pl.embedding), so it works for Visium spots, Stereo-seq bins and cellbin
alike without a tissue image. Interpretation: a cluster whose units sit together
is a spatial domain; a scattered cluster is a cell state spread across regions.
Read the top-marker maps next to the cluster map: a marker that paints a known
structure confirms the domain, one that paints nothing useful questions it.
Pitfalls: a cluster map alone can look like anatomy even when the clusters are
technical (batch, depth); confirm with markers and the tissue image when
available. This step needs the cluster step to have run.