# MAG table

**Purpose.** One row per MAG joining quality, taxonomy, abundance and
dereplication cluster — the deliverable table for reports.

**Principle.** Left joins on bin stem keep every quality-assessed bin even
when optional inputs are missing; `joined` records what actually
contributed. GTDB-Tk classification strings (`d__...;p__...;...;s__...`)
split into rank columns; empty rank names mean GTDB could not resolve that
rank, not absence of the lineage.

**Interpretation.** GTDB-Tk only classifies bins passing its
completeness/contamination gate (mag defaults 50/10) — an unclassified bin
may simply be below the gate, not novel. dRep `secondary_cluster` at ANI
0.95 approximates species; cluster representatives are the highest-quality
bin per cluster, so abundance should be summed per cluster, not per bin,
when reporting species-level profiles.

**Pitfalls.** CoverM mean coverage is depth, not abundance: divide by
genome length (or use covered fraction with a stated minimum) before
cross-genome comparison. Bins from different samples in one dRep run are
dereplicated together — the same species from two samples collapses to one
cluster; keep per-sample abundance via CoverM, not bin counts. GTDB
taxonomy names change between releases — record the database release
(r2xx) in the report.