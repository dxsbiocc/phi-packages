Purpose: quantify per-cell quality before any threshold is chosen. Principle:
n_fragment separates empty barcodes (low) from doublets/merged cells (high);
frac_dup marks over-amplified libraries; frac_mito marks damaged chromatin;
TSS enrichment (TSSe) is the ATAC-specific signal-to-noise metric — fragments
should pile up at transcription start sites with nucleosome periodicity, so a
cell with enough fragments but low TSSe is still noise. Interpretation: read the
quantiles, not a round number; TSSe needs a gene annotation (gene_anno) and is
skipped with a gap message without one. Pitfalls: a wrong mitochondrial pattern
silently matches no gene; the fragment-size distribution should show the
nucleosome-free (<100bp) and mono-nucleosome (~200bp) peaks — a flat distribution
means poor enrichment. There is no spatial coordinate here; the unit is a
barcode, not a spot.