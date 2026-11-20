Purpose: rank genes whose expression is organised in space rather than scattered
at random. Principle: spatial autocorrelation measures whether neighbouring units
(on the spatial-neighbors graph) have similar expression. Moran's I is high when
a gene paints a coherent region; Geary's C is its local-difference counterpart
(low = structured). Interpretation: read the ranked SVG table — top genes are
candidates for spatial domains or gradients, and mapping them (spatial-plots)
confirms the pattern. Pitfalls: autocorrelation depends entirely on the adjacency
graph, so a different neighbor_method/k changes the ranking; significance is a
screen, not proof of a mechanism; a gene structured only by depth or batch will
also score high, so inspect the map before believing biology.