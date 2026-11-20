Purpose: make the run reproducible and give every later cell its shared helpers.
Principle: one config cell holds every tunable value; the setup cell seeds the
RNG, fixes figure style, and defines the bookkeeping (figure_paths,
analysis_steps, the per-section result dicts) that both the notebook and the
script tools read. Interpretation: the printed scanpy/squidpy versions and seed
are the provenance of every plot and cluster below — quote them when a result
must be reproduced. Pitfalls: changing the seed changes Leiden clusters and the
UMAP; a result that will not repeat is usually a seed or version difference, not
biology. Nothing here touches the data.