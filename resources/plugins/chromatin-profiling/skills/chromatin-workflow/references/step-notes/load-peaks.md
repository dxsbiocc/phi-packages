Purpose: import called peaks per sample as a GRangesList — the pure-downstream
entry point. Principle: peak calling (MACS2 for ChIP, SEACR for CUT&Tag) already
happened upstream, in the nf-core/chipseq wrapper or the user's SEACR run; this
step only reads the resulting BED/narrowPeak/broadPeak files. narrowPeak/broadPeak
carry signal columns that rtracklayer keeps as metadata; a plain BED is three
columns. Interpretation: paths come from peak_paths, else the sample sheet's Peaks
column, else a directory of peak files, else data_path itself; relative paths
resolve against the sample sheet's directory. Each sample reports its peak count
and median width — the width is the first clue to narrow (TF) vs broad (histone)
marks. Pitfalls: this is not peak calling and never expects FASTQ; a count-table
input has no peaks, so FRiP/consensus are unavailable from it. Sample names are the
file basenames, which must match the BAM basenames for FRiP pairing.