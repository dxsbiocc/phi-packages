# ROSE super-enhancer stitching & tangent inflection (Young Lab ROSE logic in
# pure GenomicRanges + ggplot2): optionally excludes promoter peaks within TSS +/-
# tss_exclude_bp, stitches constituent peaks within stitch_distance (default
# 12,500 bp), ranks summed signal ascending, and identifies the exact tangent
# inflection point (slope = 1 on normalized [0,1]x[0,1] curve).
need("GenomicRanges", "ggplot2")

if (is.null(peak_list) || length(peak_list) == 0L) {
  stop("super-enhancer requires peak_list (run load-peaks first)", call. = FALSE)
}

peaks_gr <- peak_list[[1L]]
sample_name <- names(peak_list)[[1L]]
if (length(peaks_gr) < 3L) {
  stop(sprintf("super-enhancer needs >= 3 peaks; '%s' has %d", sample_name, length(peaks_gr)), call. = FALSE)
}

# Resolve per-constituent signal column (explicit se_signal_col -> signalValue -> score -> 1)
mcols_names <- colnames(GenomicRanges::mcols(peaks_gr))
sig_vec <- NULL
if (!is.null(se_signal_col) && nzchar(se_signal_col) && se_signal_col %in% mcols_names) {
  sig_vec <- suppressWarnings(as.numeric(GenomicRanges::mcols(peaks_gr)[[se_signal_col]]))
} else if ("signalValue" %in% mcols_names) {
  sig_vec <- suppressWarnings(as.numeric(GenomicRanges::mcols(peaks_gr)$signalValue))
} else if ("score" %in% mcols_names) {
  sig_vec <- suppressWarnings(as.numeric(GenomicRanges::mcols(peaks_gr)$score))
}
if (is.null(sig_vec) || all(is.na(sig_vec))) {
  analysis_notes <- c(analysis_notes, "super-enhancer: no numeric signalValue/score column; using peak width weighting")
  sig_vec <- as.numeric(GenomicRanges::width(peaks_gr))
}
sig_vec[is.na(sig_vec) | sig_vec < 0] <- 0
GenomicRanges::mcols(peaks_gr)$constituent_signal <- sig_vec

# Optional TSS promoter exclusion
n_excluded_tss <- 0L
if (!is.null(tss_bed) && nzchar(tss_bed) && file.exists(tss_bed)) {
  tss_gr <- rtracklayer::import(tss_bed)
  tss_win <- suppressWarnings(GenomicRanges::promoters(tss_gr, upstream = tss_exclude_bp, downstream = tss_exclude_bp))
  keep_idx <- which(GenomicRanges::countOverlaps(peaks_gr, tss_win, ignore.strand = TRUE) == 0L)
  n_excluded_tss <- as.integer(length(peaks_gr) - length(keep_idx))
  peaks_gr <- peaks_gr[keep_idx]
} else if (!is.null(txdb) || (!is.null(gtf_path) && nzchar(gtf_path) && file.exists(gtf_path))) {
  if (is.null(txdb)) {
    txdb <- GenomicFeatures::makeTxDbFromGFF(gtf_path)
  }
  genes_gr <- GenomicFeatures::genes(txdb)
  tss_win <- suppressWarnings(GenomicRanges::promoters(genes_gr, upstream = tss_exclude_bp, downstream = tss_exclude_bp))
  keep_idx <- which(GenomicRanges::countOverlaps(peaks_gr, tss_win, ignore.strand = TRUE) == 0L)
  n_excluded_tss <- as.integer(length(peaks_gr) - length(keep_idx))
  peaks_gr <- peaks_gr[keep_idx]
}

if (length(peaks_gr) < 3L) {
  stop("fewer than 3 non-TSS constituent enhancers remained after filtering", call. = FALSE)
}

# ROSE stitching: merge enhancers within stitch_distance
stitched_gr <- GenomicRanges::reduce(peaks_gr, min.gapwidth = as.integer(stitch_distance) + 1L, ignore.strand = TRUE)
hits <- GenomicRanges::findOverlaps(stitched_gr, peaks_gr, ignore.strand = TRUE)
st_idx <- S4Vectors::queryHits(hits)
pk_idx <- S4Vectors::subjectHits(hits)

n_const <- integer(length(stitched_gr))
agg_sig <- numeric(length(stitched_gr))
sig_vals <- GenomicRanges::mcols(peaks_gr)$constituent_signal
for (k in seq_along(st_idx)) {
  s_i <- st_idx[[k]]
  n_const[[s_i]] <- n_const[[s_i]] + 1L
  agg_sig[[s_i]] <- agg_sig[[s_i]] + sig_vals[[pk_idx[[k]]]]
}

ord_asc <- order(agg_sig, decreasing = FALSE)
stitched_sorted <- stitched_gr[ord_asc]
sig_sorted <- agg_sig[ord_asc]
const_sorted <- n_const[ord_asc]
n_dom <- length(sig_sorted)

max_sig <- max(sig_sorted)
if (max_sig <= 0) max_sig <- 1.0
x_norm <- seq_len(n_dom) / n_dom
y_norm <- sig_sorted / max_sig
# Exact tangent point of slope 1 on normalized [0, 1] x [0, 1] curve: maximize (x - y)
cut_idx <- which.max(x_norm - y_norm)
cutoff_signal <- sig_sorted[[cut_idx]]

is_se <- sig_sorted > cutoff_signal
if (!any(is_se) && n_dom >= 1L) is_se[[n_dom]] <- TRUE

se_df <- data.frame(
  domain_id = sprintf("%s_SE_%04d", sample_name, seq_len(n_dom)),
  chrom = as.character(GenomicRanges::seqnames(stitched_sorted)),
  start = as.integer(GenomicRanges::start(stitched_sorted)),
  end = as.integer(GenomicRanges::end(stitched_sorted)),
  width_bp = as.integer(GenomicRanges::width(stitched_sorted)),
  n_constituents = const_sorted,
  signal = round(sig_sorted, 4),
  rank_asc = seq_len(n_dom),
  rank_desc = rev(seq_len(n_dom)),
  class = ifelse(is_se, "SE", "TE"),
  stringsAsFactors = FALSE
)
se_df <- se_df[order(se_df$rank_desc), ]

out_dir <- dirname(normalizePath(se_output, mustWork = FALSE))
if (nzchar(out_dir) && !dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.table(se_df, file = se_output, sep = "\t", quote = FALSE, row.names = FALSE)

se_only_bed <- sub("[.][^.]*$", ".bed", se_output)
if (identical(se_only_bed, se_output)) se_only_bed <- paste0(se_output, ".bed")
se_sub <- se_df[se_df$class == "SE", , drop = FALSE]
bed_df <- data.frame(
  chrom = se_sub$chrom,
  start = pmax(0L, se_sub$start - 1L),
  end = se_sub$end,
  name = se_sub$domain_id,
  score = round(se_sub$signal, 2),
  stringsAsFactors = FALSE
)
utils::write.table(bed_df, file = se_only_bed, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

plot_df <- data.frame(
  rank = seq_len(n_dom),
  signal = sig_sorted,
  class = ifelse(is_se, "Super-Enhancer (SE)", "Typical Enhancer (TE)")
)
p_rose <- ggplot2::ggplot(plot_df, ggplot2::aes(x = rank, y = signal, color = class)) +
  ggplot2::geom_point(size = 1.6, alpha = 0.85) +
  ggplot2::geom_hline(yintercept = cutoff_signal, linetype = "dashed", color = "#334155") +
  ggplot2::geom_vline(xintercept = cut_idx, linetype = "dotted", color = "#64748b") +
  ggplot2::scale_color_manual(values = c("Super-Enhancer (SE)" = "#e11d48", "Typical Enhancer (TE)" = "#64748b")) +
  ggplot2::theme_minimal() +
  ggplot2::labs(
    title = sprintf("ROSE super-enhancers (%s): %d SE / %d TE", sample_name, sum(is_se), sum(!is_se)),
    x = "Enhancer rank (ascending)",
    y = "Stitched enhancer signal",
    color = NULL
  )
show_plot(p_rose, "super-enhancer-hockeystick", width = 6.8, height = 4.6)

analysis_steps <- c(analysis_steps, "super-enhancer")
se_result <- list(
  sample = sample_name,
  n_constituent_peaks = as.integer(length(peaks_gr)),
  n_excluded_tss = n_excluded_tss,
  n_stitched_domains = as.integer(n_dom),
  n_super_enhancers = as.integer(sum(is_se)),
  n_typical_enhancers = as.integer(sum(!is_se)),
  stitch_distance_bp = as.integer(stitch_distance),
  signal_cutoff = round(cutoff_signal, 4),
  output_table = se_output,
  output_bed = se_only_bed,
  top_super_enhancers = utils::head(se_sub, 15L)
)
cat(sprintf("super-enhancer (%s): %d SEs from %d stitched domains (cutoff=%.2f)\n",
            sample_name, sum(is_se), n_dom, cutoff_signal))
