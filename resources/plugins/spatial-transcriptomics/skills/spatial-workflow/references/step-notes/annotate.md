Purpose: score each domain against a marker panel to propose what it is.
Principle: score = panel mean minus the mean of everything else, so panel sizes
stay comparable; this is a heuristic ranking, not a trained classifier.
Interpretation: a label needs positive markers AND the absence of competing-type
markers — read the dot plot before accepting the automatic call. The built-in
panel is a generic example; spatial regions often need tissue-specific markers,
so pass panel_path with a curated JSON. Pitfalls: a spatial spot/bin mixes
several cells, so a high score means "this region is enriched for", not "this
unit is one cell of that type"; absent panel genes (human MT- vs mouse mt-)
silently weaken a score; a domain with no evidence should stay unassigned rather
than be forced into a label.