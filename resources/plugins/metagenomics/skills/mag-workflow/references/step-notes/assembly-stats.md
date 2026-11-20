# Assembly statistics

**Purpose.** Decide whether an assembly is worth binning, and compare
assembler/co-assembly choices.

**Principle.** N50 = contig length such that half the assembly lies in
contigs ≥ that length; L50 = how many contigs reach it. QUAST reports these
at several minimum-length cutoffs (>=0, >=1000, >=5000 bp); the tool reads
the >=0 bp rows so nothing is hidden by cutoff choice.

**Interpretation.** Metagenome assemblies are mixtures of genomes at
uneven depth — N50 of 2–10 kb is ordinary for complex communities, and
total length reflects the community's genome sum plus strain variation.
MEGAHIT vs metaSPAdes vs co-assembly differences of 2× in N50 are common
and not automatically "better/worse"; binning yield is the real test.

**Pitfalls.** Do not apply isolate expectations (N50 > 100 kb) to
metagenomes. GC% far from the expected community range hints at
contamination or organelle/host sequence. QUAST without a reference gives
no misassembly count — fragmentation stats alone cannot certify correctness.