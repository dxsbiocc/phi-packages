Purpose: persist the analysis so a later session or script resumes from here.
Principle: the RDS holds the DBA object (or the consensus GRanges when there is no
DBA); the tables beside it are the reviewable, language-neutral outputs.
Interpretation: it writes consensus.bed (rtracklayer::export), counts.csv,
normalized_counts.csv, occupancy.csv, diff_binding.csv (DiffBind) or
csaw_regions.csv (csaw), and annotation.csv into the output directory, creating it
if needed. A saved DBA can be re-inspected (chip_inspect reports saved-object) and
resume at diff/annotate. Pitfalls: with no DBA object and no consensus peakset there
is nothing to save as RDS (a note says so) — tables may still be written. The RDS
is R-specific; the CSV/BED tables are the portable deliverables to quote in the
answer. Saving does not close any backed object (R has none here, unlike the
snapatac2 scATAC route).