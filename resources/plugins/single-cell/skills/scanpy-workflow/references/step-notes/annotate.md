Purpose: score clusters against marker panels so labels rest on stated evidence.
Principle: score = panel mean minus the mean of everything else, keeping panel
sizes comparable (the same rule as the Seurat route); a real identification
needs positive markers AND the absence of competing types' markers.
Interpretation: read the score margins and the dot plot — a narrow win between
two similar panels is not an identification; work coarse-to-fine (lineage first,
subtype second). Pitfalls: a reference panel forces a label onto every cluster
even when the type is absent, so an unscored cluster stays unassigned; zero
usable panel genes means an identifier mismatch (human MT- vs mouse mt-,
Ensembl IDs), not a biology result; CellTypist needs log1p 1e4 input and a
project environment.
