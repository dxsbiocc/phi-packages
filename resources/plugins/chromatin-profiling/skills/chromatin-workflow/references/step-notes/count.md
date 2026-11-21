Purpose: count reads in the consensus peaks for every sample — the affinity matrix
the differential analysis tests. Principle: GenomicAlignments::summarizeOverlaps
counts fragments overlapping each peak in Union mode (a read counts once even if
it overlaps two peaks), with inter.feature=FALSE so shared reads are not discarded.
When a DBA object exists, DiffBind::dba.count(score=DBA_SCORE_READS) also populates
it so dba.analyze can run. Interpretation: it reports the matrix shape and
per-sample library sizes; a vastly smaller library is under-sequenced, a vastly
larger one may be over-amplified. Pitfalls: needs indexed BAM and a consensus
peakset; paired_end must match the library (single-end counts reads, paired-end
counts fragments) or every peak is miscounted. These are raw counts — the
differential engines normalize internally; the normalized copy from the next step
is for exploration only.