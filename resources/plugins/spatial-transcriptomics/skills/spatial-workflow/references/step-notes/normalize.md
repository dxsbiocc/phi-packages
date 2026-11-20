Purpose: make counts comparable across spatial units. Principle: library size
differs by orders of magnitude between spots/bins, so normalize_total rescales
each unit to scale_factor counts and log1p stabilizes variance so a few abundant
genes do not dominate every distance; the pre-scale result is snapshotted to the
lognorm layer. Interpretation: counts keeps the raw UMIs (for HVG and
deconvolution), lognorm keeps the expression of record (for markers, panel
scoring, SVG and ligrec). Pitfalls: normalizing an already-normalized object
double-transforms it — re-import raw counts instead; scaled values (next steps)
are z-scores and must never enter a marker test or a fold change.