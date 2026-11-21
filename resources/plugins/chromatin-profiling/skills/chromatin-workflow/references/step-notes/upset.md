Purpose: peak overlap across samples as an UpSet plot — how many sites are private
to one sample versus shared by two, three, or all. Principle:
ComplexHeatmap::make_comb_mat builds the combination matrix from the consensus
occupancy matrix, and UpSet draws the intersection sizes; it is the visual
companion to the reproducibility overlap fractions. Interpretation: a condition
whose replicates share most peaks is reproducible; samples that share almost
nothing differ biologically or technically. Pitfalls: it needs a consensus
occupancy matrix with >=2 samples; failing that it falls back to
ChIPseeker::upsetplot on a single annotated csAnno when ggupset is installed (a
ChIPseeker Suggests, not shipped), or skips with a note. A count-table
input has no occupancy matrix, so the fallback or skip applies — the absence is
not a clean result. UpSet shows set structure, not significance.