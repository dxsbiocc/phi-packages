Purpose: test which cell types/domains sit next to each other more than chance.
Principle: nhood_enrichment compares the observed count of neighbour pairs
against a permutation null, giving a z-score per (cluster_i, cluster_j);
co_occurrence measures how often two types appear within a spatial interval,
revealing gradients. Interpretation: a positive z-score means two domains are
adjacent more than expected (an interface, a niche); read it together with the
spatial map. Pitfalls: enrichment is about adjacency, not interaction — neighbours
need not signal to each other; the result inherits the adjacency graph and the
cluster labels, so noisy clusters or a wrong graph produce spurious niches; with
few units the permutation null is weak.