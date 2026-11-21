# Annotate peaks to genomic features with ChIPseeker: promoter / exon / intron /
# intergenic distribution, nearest gene, and distance to the TSS. The lean
# annotation strategy means the target-genome TxDb comes from the user — a GTF /
# GFF3 (makeTxDbFromGFF) or a prebuilt TxDb .sqlite (loadDb). The plugin bundles no
# target-genome TxDb; without one this step reports a gap instead of guessing.
if (!is.null(txdb_path) && nzchar(txdb_path)) {
  if (!file.exists(txdb_path)) stop("TxDb file not found: ", txdb_path, call. = FALSE)
  txdb <- AnnotationDbi::loadDb(txdb_path)
  cat("loaded user TxDb:", txdb_path, "\n")
} else if (!is.null(gtf_path) && nzchar(gtf_path)) {
  if (!file.exists(gtf_path)) stop("annotation file not found: ", gtf_path, call. = FALSE)
  cat("building TxDb from", gtf_path, "(this can take a minute for a full GTF)\n")
  txdb <- txdbmaker::makeTxDbFromGFF(gtf_path)
} else {
  stop("annotation needs a user GTF/GFF3 (gtf_path) or TxDb .sqlite (txdb_path). ",
       "The plugin bundles no target-genome TxDb (lean annotation strategy).", call. = FALSE)
}

peaks <- if (!is.null(consensus_peaks)) {
  consensus_peaks
} else if (!is.null(peak_list)) {
  as.list(peak_list)
} else {
  stop("load peaks or build a consensus before annotating", call. = FALSE)
}

anno_args <- list(TxDb = txdb, tssRegion = c(-tss_upstream, tss_downstream),
                  level = "gene", verbose = FALSE)
# Gene symbols need a user org.db; without it annotation still gives gene ids,
# features and TSS distance. The bundled hg19 TxDb (a ChIPseeker dependency) is
# never used as the annotation source for the target genome.
if (!is.null(org_db) && nzchar(org_db)) {
  if (have(org_db)) {
    anno_args$annoDb <- org_db
  } else {
    cat("note: ", org_db, " is not installed; gene symbols unavailable, reporting gene ids. ",
        "A project environment (env_request) can provide it.\n", sep = "")
  }
}
# annotatePeak takes one peak set (a GRanges or a file), not a list: annotate each
# sample separately so plotAnnoBar / plotDistToTSS can compare them side by side.
annotation_result <- if (inherits(peaks, "GRanges")) {
  do.call(ChIPseeker::annotatePeak, c(list(peak = peaks), anno_args))
} else {
  lapply(peaks, function(p) do.call(ChIPseeker::annotatePeak, c(list(peak = p), anno_args)))
}

show_plot(ChIPseeker::plotAnnoBar(annotation_result), "feature-distribution")
show_plot(ChIPseeker::plotDistToTSS(annotation_result), "tss-distance")

anno_df <- if (inherits(annotation_result, "csAnno")) {
  as.data.frame(annotation_result)
} else {
  do.call(rbind, lapply(names(annotation_result), function(nm) {
    df <- as.data.frame(annotation_result[[nm]])
    df$sample <- nm
    df
  }))
}
cat(sprintf("annotated %d peaks; feature columns: %s\n",
            nrow(anno_df), paste(head(colnames(anno_df), 8), collapse = ", ")))
analysis_steps <- c(analysis_steps, "annotate")