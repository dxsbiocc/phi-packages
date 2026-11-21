Purpose: annotate peaks to genomic features with ChIPseeker — promoter/exon/intron/
intergenic distribution, nearest gene, and distance to the TSS. Principle: the lean
annotation strategy means the target-genome TxDb comes from the user — a GTF/GFF3
(txdbmaker::makeTxDbFromGFF) or a prebuilt TxDb .sqlite (AnnotationDbi::loadDb); the
plugin bundles no target-genome TxDb, and without one the step reports a gap instead
of guessing.
Interpretation: annotatePeak(peak, TxDb, tssRegion=c(-up,down), level="gene")
returns a csAnno for one peak set; with per-sample peaks the fragment lapplys it
into a named list of csAnno. feature_distribution
tabulates the annotation column. A TF set enriched at promoters suggests
transcriptional regulation; a broad mark over gene bodies suggests elongation/
heterochromatin. Pitfalls: the hg19 TxDb that ships as a ChIPseeker conda dependency
is a build-time artifact, NEVER the annotation source for hg38/mm10 peaks. Gene
symbols need a user org.db (annoDb); without it annotation gives gene ids, not
symbols (gene_symbols=false). Annotation is positional inference, not measured
target-gene regulation; motif enrichment (HOMER/MEME) is a separate gap.