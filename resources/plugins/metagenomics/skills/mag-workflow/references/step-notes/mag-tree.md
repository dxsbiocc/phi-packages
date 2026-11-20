# MAG tree distances

**Purpose.** Quantify where MAGs sit on the GTDB reference tree: pairwise
patristic distances between user MAGs, and each MAG's nearest GTDB reference
tip with its distance — without a phylogenetics library.

**Principle.** GTDB-Tk `classify.tree` is a Newick tree placing user genomes
into the GTDB backbone (tens of thousands of reference tips). Patristic
distance = sum of branch lengths along the path between two tips, a
marker-gene divergence measure. The parser is pure Python (arrays, no
Bio.Phylo); user tips come from `--summary`'s user_genome column when given
(exact, recommended), else a heuristic excluding `RS_/GB_ GCF_/GCA_`
reference names — the JSON reports which was used. Nearest-reference search
walks each reference tip up to the root once and takes the minimum over
shared ancestors, so cost scales with reference count, not MAG × reference
pairs.

**Interpretation.** A small distance to a named reference supports the
GTDB-Tk species call; a large distance to every reference is consistent
with a novel lineage — but check the completeness/contamination gate first
(incomplete bins place unpredictably). The distance matrix clusters MAGs
by marker-gene divergence; it is not an ANI matrix and does not replace
dRep for dereplication decisions.

**Pitfalls.** Patristic ≠ ANI: branch lengths are substitutions per site on
the backbone alignment, not genome-wide identity. No tree figure is
rendered here (no phylogenetics library in the environment); the distance
matrix is the deliverable. Trees without user tips fail loudly — pass the
matching `--summary`.