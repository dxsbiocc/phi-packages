Purpose: persist the object so the next step or a colleague can continue without
recomputing. Principle: the AnnData is written to .h5ad with its counts and lognorm
layers, spatial coordinates in obsm, the adjacency graph in obsp, and all obs
annotations intact; the JSON report records the parameters used. Interpretation:
the saved object is the audit record — reload it and the shape, layers and
annotations match what the report claimed. Pitfalls: save only what you intend to
carry forward (a scaled X is not the layer of record for markers); the output
path must be writable and is recorded relative to the working directory; an
object saved without coordinates cannot later run spatial statistics.