Purpose: resolve the blacklist of artifactual high-signal regions to exclude
from feature selection and peak calling. Principle: ENCODE blacklist regions
(satellite repeats, gaps, anomalous pileups) attract fragments regardless of
biology, so leaving them in creates spurious "accessible" features and peaks.
Interpretation: the shipped ENCODE v2 blacklists for hg38 (636 regions) and mm10
(3435 regions) are inlined offline as base64-gzip and decoded to a temp BED; a
user blacklist path overrides them. For a genome with no shipped blacklist
(e.g. hg19) the step returns none and says so. Pitfalls: a blacklist is
genome-build specific — hg38 regions are wrong on mm10; passing the wrong build
silently removes the wrong loci. Blacklisting is not optional cleanup: it changes
which features survive selection and which peaks MACS3 calls. Source: Amemiya,
Kundaje & Boyle, Sci Rep 9:9354 (2019), doi:10.1038/s41598-019-45839-z.