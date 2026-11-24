Purpose: partition tissue into recurring microenvironments (spatial niches) by
clustering units that share similar spatial neighbourhoods, and optionally track
niche proportion shifts across conditions or disease stages. Principle: builds a
row-normalized spatial aggregation operator (A + I)^hops on the spatial graph
(a pure-conda UTAG-style message-passing step), averages either local cell-type
composition (or deconvolution proportions) or lognorm expression, and clusters
the aggregated representation with Leiden. Interpretation: read the spatial
niche map alongside the niche composition heatmap to name each tissue domain
(e.g. tumour core, invasive margin, lymphoid follicle), and check per-condition
niche shifts when group_column is set. Pitfalls: niche boundaries depend on the
spatial graph, niche_hops, and niche_resolution; composition mode inherits the
accuracy of group_by or deconvolution; cross-condition shifts are descriptive
unless supported by biological replicates.
