Purpose: test which cell types/domains sit next to each other more than chance,
quantify domain contact frequencies, spatial graph centrality, and point-pattern
clustering across scales. Principle: nhood_enrichment compares observed neighbour
pairs against a permutation null (z-score); co_occurrence measures probability
within growing spatial radii; interaction_matrix reports normalized contact
proportions; centrality_scores computes degree, clustering, and closeness
centrality on the spatial graph; ripley (L/F/G) tests whether each domain is
spatially clustered or dispersed across distances. Interpretation: combine
enrichment z-scores, contact proportions, centrality (core vs dispersed/boundary
domains), and Ripley's L(r) - r (> 0 clustered, < 0 dispersed) with the spatial
map. Pitfalls: adjacency is not signalling; all five metrics inherit the spatial
graph and cluster labels, so noisy clusters or a poorly chosen graph produce
spurious patterns; with very few units permutations and Ripley envelopes are weak.
