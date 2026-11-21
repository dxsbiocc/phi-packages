# Count reads in the consensus peaks for every sample — the affinity matrix the
# differential analysis tests. summarizeOverlaps counts fragments overlapping each
# peak (Union mode: a read counts once even if it overlaps two peaks). When a DBA
# object exists, DiffBind's dba.count also populates it so dba.analyze can run.
if (is.null(consensus_peaks)) stop("build a consensus peakset first (consensus)", call. = FALSE)
if (is.null(bam_files) || length(bam_files) == 0) stop("load BAM first (load-bam)", call. = FALSE)

single_end <- !isTRUE(paired_end)
overlaps <- GenomicAlignments::summarizeOverlaps(
  consensus_peaks, bam_files, mode = "Union",
  singleEnd = single_end, ignore.strand = TRUE, inter.feature = FALSE)
count_matrix <- SummarizedExperiment::assay(overlaps)
rownames(count_matrix) <- paste0(seqnames(consensus_peaks), ":",
                                 start(consensus_peaks), "-", end(consensus_peaks))
colnames(count_matrix) <- tools::file_path_sans_ext(basename(bam_paths))

if (!is.null(dba_object)) {
  dba_object <- DiffBind::dba.count(dba_object, peaks = consensus_peaks,
                                  score = DiffBind::DBA_SCORE_READS)
  cat("DBA object counted over the consensus peakset\n")
}

cat(sprintf("count matrix: %d peaks x %d samples\n", nrow(count_matrix), ncol(count_matrix)))
cat("library sizes:", paste(format(colSums(count_matrix), big.mark = ","), collapse = ", "), "\n")
analysis_steps <- c(analysis_steps, "count")