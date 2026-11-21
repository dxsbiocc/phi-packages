Purpose: ask which TF motifs are over-represented in each cluster's marker
regions, turning open chromatin into a regulatory hypothesis. Principle:
motif_enrichment scans the marker region sequences (extracted from a genome
FASTA) against a motif database and tests over-representation with a binomial or
hypergeometric test against a background. Interpretation: an enriched motif
suggests a TF family active in that cluster's regulatory landscape; read the
per-cluster enrichment tables alongside the marker genes. Pitfalls: this is
sequence-level PREDICTION, not measured TF binding — motifs are short and
degenerate, many TFs share motifs, and a hit needs experimental validation; it
needs a genome FASTA (genome_fasta) and a motif database (motif_db, or snapatac2
cis_bp which downloads), without genome_fasta the step reports a gap; the
background and database choice change the result, so state them. chromVAR
per-cell TF deviations are a different, R-route analysis NOT implemented here —
do not call motif enrichment "chromVAR".