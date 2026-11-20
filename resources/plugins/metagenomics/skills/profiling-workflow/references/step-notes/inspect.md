# Inspect

**Purpose.** Identify what a file actually is before any analysis.

**Principle.** Detection is header- and shape-based: Kraken2 reports have six
tab columns with a rank code in column 4; Bracken adds `new_est_reads` /
`taxonomy_lvl`; MetaPhlAn profiles carry `clade_name` + `relative_abundance`
with `k__|p__|...|s__` lineages; HUMAnN tables start with `# Gene Family` or
pathway IDs and may be stratified (`gene|bug`).

**Interpretation.** `ready_steps` lists what the input supports. `sample`
comes from the file stem — rename files to real sample names before merging
if stems are lane IDs.

**Pitfalls.** Kraken2 sub-ranks (S1, G2) are not the plain rank code; merge
keeps exact matches only. A "matrix" with non-numeric columns is rejected
downstream — fix the file, do not coerce silently.