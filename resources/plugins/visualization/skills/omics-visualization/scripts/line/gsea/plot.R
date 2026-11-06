#!/usr/bin/env Rscript

# Template-ID: line-gsea
#
# Purpose:
#   GSEA running-enrichment plot (enrichplot::gseaplot2 style): one
#   running-score curve per gene set, a hit-barcode row per set, and
#   the ranked list metric as a grey area, stacked on a shared rank
#   axis. A single-set classic view adds a red-blue rank strip, the
#   leading-edge line at the peak and supplied NES / FDR.
#
# Inputs:
#   ranked.tsv   - gene, metric (full ranked list; sorted descending
#                  by the script, ties keep file order)
#   Sidecars beside the input (filenames in CONFIG):
#     genesets.tsv  - set, gene
#     running.tsv   - set, rank, score (supplied curve; e.g. from
#                     enrichplot:::gsInfo or fgsea plotEnrichment data)
#     set_stats.tsv - optional set, NES, padj (labels only)
#
# Output:
#   A PDF, PNG, or SVG GSEA enrichment plot.
#
# Dependencies:
#   ggplot2, readr, patchwork
#
# Example:
#   Rscript plot.R ranked.tsv output.pdf
#
# Agent adaptation:
#   Edit only CONFIG (sets, running source, palette, panel heights,
#   metric panel, rank strip, leading edge, stats labels). Extra sets
#   are extra rows; the single-set classic view is a variant, not a
#   new id. Do not run GSEA, permutations, NES, p or FDR here.
#
# Scientific assumptions:
#   The ranked metric, gene sets, running scores and NES / FDR are
#   supplied. running = "compute" redraws the deterministic weighted
#   running sum (Subramanian 2005, weight p) from the ranked list for
#   display only; it does not test or normalise anything.

local({
    file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    if (!length(file_arg)) {
        stop("Run this file with Rscript.", call. = FALSE)
    }
    dir <- dirname(normalizePath(sub("^--file=", "", file_arg[[1]])))
    for (i in seq_len(8)) {
        candidate <- file.path(dir, "lib", "common.R")
        if (file.exists(candidate)) {
            source(normalizePath(candidate), chdir = FALSE)
            return(invisible())
        }
        parent <- dirname(dir)
        if (identical(parent, dir)) break
        dir <- parent
    }
    stop("Cannot find scripts/lib/common.R", call. = FALSE)
})

io <- parse_io_args()

# -----------------------------------------------------------------------------
# CONFIG  (edit this block for a new dataset)
# -----------------------------------------------------------------------------
config <- list(
    columns = list(gene = "gene", metric = "metric"),
    sidecars = list(
        genesets = "genesets.tsv",
        running = "running.tsv",
        stats = "set_stats.tsv"
    ),
    sets = NULL,                 # NULL = all sets in genesets.tsv order
    running = "supplied",        # "supplied" (running.tsv) or "compute"
    weight = 1,                  # running-sum weight when computing
    palette = "Qualitative.Bold",   # or "Family.Name:index" for one colour
    heights = c(curve = 3, barcode = 0.28, metric = 1.6),  # barcode per set
    show_metric = TRUE,
    rank_strip = FALSE,          # red-blue rank colour strip under barcodes
    leading_edge = FALSE,        # dashed line at each curve's peak
    annotate_stats = FALSE,      # NES / FDR (legend, or corner if one set)
    legend_inside = c(0.99, 0.99),
    labels = list(
        x = "Rank in Ordered Dataset",
        curve = "Running Enrichment Score",
        metric = "Ranked List Metric"
    ),
    size = list(width = 9.5, height = 6.5)
)

load_packages(c("ggplot2", "readr", "patchwork"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
side <- function(key, cols, required = TRUE) {
    fn <- config$sidecars[[key]]
    path <- if (is.null(fn)) NA else file.path(dirname(io$input), fn)
    if (is.na(path) || !file.exists(path)) {
        if (required) stop(paste("Missing sidecar:", fn), call. = FALSE)
        return(NULL)
    }
    d <- as.data.frame(read_table_auto(path), stringsAsFactors = FALSE)
    require_columns(d, as.list(stats::setNames(cols, cols)))
    d
}

rk <- as.data.frame(read_table_auto(io$input), stringsAsFactors = FALSE)
require_columns(rk, config$columns)
rk <- data.frame(gene = as.character(rk[[config$columns$gene]]),
    metric = as.numeric(rk[[config$columns$metric]]),
    stringsAsFactors = FALSE)
if (any(!is.finite(rk$metric))) stop("metric must be finite.", call. = FALSE)
if (anyDuplicated(rk$gene)) stop("gene must be unique.", call. = FALSE)
rk <- rk[order(-rk$metric), , drop = FALSE]
rk$rank <- seq_len(nrow(rk))
N <- nrow(rk)

gs <- side("genesets", c("set", "gene"))
sets <- if (is.null(config$sets)) unique(gs$set) else config$sets
if (!all(sets %in% gs$set)) {
    stop("CONFIG sets missing from genesets.tsv.", call. = FALSE)
}
hits <- lapply(sets, function(s) {
    r <- rk$rank[match(gs$gene[gs$set == s], rk$gene)]
    sort(r[!is.na(r)])
})
names(hits) <- sets

if (identical(config$running, "compute")) {
    run <- do.call(rbind, lapply(sets, function(s) {
        h <- hits[[s]]
        w <- abs(rk$metric[h])^config$weight
        step <- rep(-1 / (N - length(h)), N)
        step[h] <- w / sum(w)
        sc <- cumsum(step)
        data.frame(set = s, rank = c(0, seq_len(N)), score = c(0, sc))
    }))
} else {
    run <- side("running", c("set", "rank", "score"))
    run$rank <- as.numeric(run$rank)
    run$score <- as.numeric(run$score)
    run <- run[run$set %in% sets, , drop = FALSE]
    if (!all(sets %in% run$set)) {
        stop("running.tsv lacks a selected set.", call. = FALSE)
    }
}

st <- side("stats", c("set", "NES", "padj"), required = FALSE)
lab <- stats::setNames(sets, sets)
if (isTRUE(config$annotate_stats) && !is.null(st) && length(sets) > 1L) {
    i <- match(sets, st$set)
    lab <- stats::setNames(ifelse(is.na(i), sets, sprintf(
        "%s  (NES %.2f, FDR %s)", sets, st$NES[i],
        formatC(st$padj[i], format = "e", digits = 1))), sets)
}
pal <- strsplit(config$palette, ":", fixed = TRUE)[[1]]
cols <- if (length(pal) > 1L) {
    rep(palette_colors(pal[[1]])[[as.integer(pal[[2]])]], length(sets))
} else {
    palette_colors(pal[[1]], n = length(sets))
}
names(cols) <- sets
run$set <- factor(run$set, levels = sets)
peak <- do.call(rbind, lapply(sets, function(s) {
    d <- run[run$set == s, , drop = FALSE]
    d[which.max(abs(d$score)), c("set", "rank", "score")]
}))

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
xs <- list(scale_x_continuous(expand = c(0, 0)),
    coord_cartesian(xlim = c(0, N)))
th <- theme_bw(base_size = 11) +
    theme(
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.text = element_text(colour = "black"),
        plot.background = element_blank(),
        plot.margin = margin(0, 6, 0, 6)
    )
no_x <- theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    axis.title.x = element_blank())

p1 <- ggplot(run, aes(x = rank, y = score, colour = set)) +
    geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
    geom_line(linewidth = 0.9) +
    scale_colour_manual(values = cols, labels = lab, name = NULL) +
    xs + labs(y = config$labels$curve) + th + no_x +
    theme(legend.position = "inside",
        legend.position.inside = config$legend_inside,
        legend.justification = config$legend_inside,
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(size = 8))
if (length(sets) == 1L) p1 <- p1 + theme(legend.position = "none")
if (isTRUE(config$leading_edge)) {
    p1 <- p1 + geom_vline(data = peak, aes(xintercept = rank, colour = set),
        linetype = "dashed", linewidth = 0.5, show.legend = FALSE)
}
if (isTRUE(config$annotate_stats) && !is.null(st) && length(sets) == 1L) {
    i <- match(sets, st$set)
    if (!is.na(i)) {
        p1 <- p1 + annotate("text", x = N * 0.02, y = -Inf, vjust = -0.6,
            hjust = 0, size = 3.8, label = sprintf("NES = %.2f\nFDR = %s",
                st$NES[i], formatC(st$padj[i], format = "e", digits = 1)))
    }
}

bc <- do.call(rbind, lapply(seq_along(sets), function(k) {
    data.frame(set = sets[[k]], rank = hits[[k]],
        row = length(sets) - k + 1)
}))
bc$set <- factor(bc$set, levels = sets)
bar_col <- if (length(sets) == 1L) {
    stats::setNames("black", sets)
} else {
    cols
}
p2 <- ggplot(bc) +
    geom_segment(aes(x = rank, xend = rank, y = row - 0.42, yend = row + 0.42,
        colour = set), linewidth = 0.35, show.legend = FALSE) +
    scale_colour_manual(values = bar_col)
ylo <- 0.5
if (isTRUE(config$rank_strip)) {
    nb <- 200
    br <- cut(rk$rank, nb, labels = FALSE)
    strip <- data.frame(x = tapply(rk$rank, br, mean),
        m = tapply(rk$metric, br, mean))
    lim <- max(abs(strip$m))
    p2 <- p2 + geom_tile(data = strip, aes(x = x, y = 0, fill = m),
        width = N / nb, height = 0.7) +
        scale_fill_gradientn(colours = rev(palette_colors("Diverging.RdBu")),
            limits = c(-lim, lim), guide = "none")
    ylo <- -0.4
}
if (isTRUE(config$leading_edge)) {
    p2 <- p2 + geom_vline(data = peak, aes(xintercept = rank),
        colour = "grey40", linetype = "dashed", linewidth = 0.4)
}
p2 <- p2 + xs +
    scale_y_continuous(limits = c(ylo, length(sets) + 0.5), expand = c(0, 0)) +
    th + no_x + theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), axis.title.y = element_blank(),
        panel.grid = element_blank())

h <- config$heights
n_bar <- max(length(sets), 1) * h[["barcode"]] +
    if (isTRUE(config$rank_strip)) h[["barcode"]] else 0
if (isTRUE(config$show_metric)) {
    p3 <- ggplot(rk, aes(x = rank, y = metric)) +
        geom_area(fill = "grey65") +
        geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) +
        xs + labs(x = config$labels$x, y = config$labels$metric) + th
    p <- p1 / p2 / p3 + patchwork::plot_layout(
        heights = c(h[["curve"]], n_bar, h[["metric"]]))
} else {
    p2 <- p2 + theme(axis.text.x = element_text(colour = "black"),
        axis.ticks.x = element_line()) + labs(x = config$labels$x)
    p <- p1 / p2 + patchwork::plot_layout(heights = c(h[["curve"]], n_bar))
}

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output, width = config$size$width,
    height = config$size$height)
