# Coverage profile around the TSS and a peak heatmap: the shape of the signal,
# not just its location. A sharp TF peak concentrates at the motif; a broad mark
# spreads over the gene body. ChIPseeker builds a tag matrix from the peaks and
# the user TxDb, then averages it (plotAvgProf2) and draws per-peak rows
# (peakHeatmap). bigWig tracks, when present, give the same view genome-wide.
if (is.null(txdb)) {
  if (!is.null(txdb_path) && nzchar(txdb_path)) {
    txdb <- AnnotationDbi::loadDb(txdb_path)
  } else if (!is.null(gtf_path) && nzchar(gtf_path)) {
    txdb <- txdbmaker::makeTxDbFromGFF(gtf_path)
  } else {
    stop("coverage profiles need a user TxDb (gtf_path or txdb_path); run annotate first",
         call. = FALSE)
  }
}
peaks <- if (!is.null(consensus_peaks)) {
  consensus_peaks
} else if (!is.null(peak_list)) {
  as.list(peak_list)
} else {
  stop("load peaks before drawing coverage profiles", call. = FALSE)
}

profile <- ChIPseeker::plotAvgProf2(peaks, TxDb = txdb,
                                    upstream = tss_upstream, downstream = tss_downstream,
                                    conf = 0.95, res = 500,
                                    xlab = "Distance to TSS (bp)", ylab = "Read count frequency")
show_plot(profile, "tss-profile")
heatmap <- ChIPseeker::peakHeatmap(peaks, TxDb = txdb,
                                   upstream = tss_upstream, downstream = tss_downstream)
show_plot(heatmap, "peak-heatmap")
cat(sprintf("TSS profile + heatmap over %d bp upstream / %d bp downstream\n",
            tss_upstream, tss_downstream))
analysis_steps <- c(analysis_steps, "coverage-profile")