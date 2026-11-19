Purpose: flag two-cells-in-one-droplet artifacts before they invent clusters.
Principle: Scrublet simulates artificial doublets from the observed counts and
scores each cell against them; heterotypic doublets (two different cell types)
masquerade as novel "hybrid" states, homotypic ones are nearly invisible.
Interpretation: the expected rate scales with loaded cells (~4% at ~9k loaded);
inspect the scores on the UMAP and run detection per sample when libraries were
captured separately. Pitfalls: the filter cell removes flagged cells — scoring
alone changes nothing; the rate is a loading estimate, not a free knob, and
 Scrublet can fail on tiny or sparse inputs, which the report states explicitly.
