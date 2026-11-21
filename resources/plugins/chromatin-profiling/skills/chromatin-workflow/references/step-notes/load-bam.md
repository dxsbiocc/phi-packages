Purpose: open the aligned-read BAM files used for FRiP, read counting and coverage
profiles. Principle: alignment (bowtie2/BWA) already happened upstream in the
nf-core/chipseq wrapper; this step only opens the BAMs via Rsamtools::BamFileList.
Counting and FRiP need a .bai index beside each BAM. Interpretation: paths come
from bam_paths, else the sample sheet's bamReads column, else bam_dir/data_path;
relative paths resolve against the sheet directory. A missing index produces a
warning naming the files, not a silent skip. Pitfalls: samtools index is upstream
work — the plugin does not build indexes; without .bai, summarizeOverlaps/countBam
cannot run, so FRiP and counting are blocked. paired_end must match the library or
fragments are miscounted downstream.