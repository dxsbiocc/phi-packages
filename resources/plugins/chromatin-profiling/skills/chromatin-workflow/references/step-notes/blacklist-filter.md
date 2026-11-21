Purpose: remove artifactual high-signal regions before consensus and counting.
Principle: ENCODE blacklist regions (satellite repeats, gaps, anomalous pileups)
attract reads regardless of biology, so leaving them in creates spurious consensus
peaks and inflates the multiple-testing universe. Interpretation: the shipped
ENCODE v2 blacklists for hg38 (636 regions) and mm10 (3435 regions) are inlined
offline as base64-gzip and decoded to a temp BED (memDecompress + jsonlite); a user
blacklist path overrides them, and `genome` selects the build. With neither a
genome nor a blacklist the step skips and says so. blacklist_removed records how
many peaks were dropped per sample and from the consensus. Pitfalls: a blacklist
is genome-build specific — hg38 regions are wrong on mm10; passing the wrong build
silently removes the wrong loci. Blacklisting is not optional cleanup: it changes
which peaks survive into the consensus and the differential test. Source: Amemiya,
Kundaje & Boyle, Sci Rep 9:9354 (2019), doi:10.1038/s41598-019-45839-z.