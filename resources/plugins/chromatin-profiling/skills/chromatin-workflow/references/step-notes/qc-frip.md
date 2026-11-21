Purpose: compute FRiP (fraction of reads in peaks), the headline ChIP/CUT&Tag
quality metric. Principle: a low FRiP means the antibody did not enrich — most
reads are background. It is computed from BAM + peaks via Rsamtools::countBam
(mapped reads) and GenomicAlignments::summarizeOverlaps (reads in peaks, Union
mode), so it is independent of any peak caller's own report. Interpretation:
ENCODE guidance — TF ChIP ideally >1%, broad marks lower, CUT&Tag often much higher
because background is minimal. Read the per-sample frip distribution; a sample far
below its peers is a failed enrichment, not a biology difference. Pitfalls: frip
is NA with a "no matching peak set" note when BAM and peak names do not align
(positional fallback is tried first) — fix the pairing, do not read NA as 0. Needs
indexed BAM and loaded peaks; paired_end changes how fragments are counted.
Source: Landt et al., Genome Res 22:1813 (2012), doi:10.1101/gr.136184.111.