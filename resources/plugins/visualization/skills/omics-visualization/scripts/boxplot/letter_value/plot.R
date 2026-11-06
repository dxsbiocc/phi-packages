#!/usr/bin/env Rscript

# Template-ID: boxplot-letter-value
#
# Purpose:
#   Letter-value (boxen) plot: per category, nested boxes at the
#   quartiles, octiles, sixteenths, ... of the supplied values, each
#   outer level narrower and lighter, a median line, and diamond
#   outliers beyond the outermost level. Larger groups get more
#   levels, so tails stay visible without a violin.
#
# Inputs:
#   One row per observation. Default example:
#     - subtype: category (label shown as written, e.g. "Calb2*")
#     - depth: numeric value
#
# Output:
#   A PDF, PNG, or SVG letter-value plot.
#
# Dependencies:
#   ggplot2, readr
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   Edit only CONFIG (columns, depth rule or fixed level count,
#   palette, reversed y, category order, outlier shape, labels).
#   Extra categories are extra rows. Reversed y and order are CONFIG,
#   not new ids. A plain quartile box is boxplot-basic; a density
#   outline is boxplot-violin.
#
# Scientific assumptions:
#   Values are supplied. Letter values are display quantiles of the
#   supplied values (as boxplot quartiles are); no test is run.
#   Level count k = floor(log2(n)) - 3 (Tukey rule), at least 1.

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
    columns = list(x = "subtype", y = "depth"),
    depth = "tukey",             # "tukey" or an integer level count
    max_levels = 8,
    box_width = 0.8,
    min_width = 0.18,            # outermost width, fraction of box_width
    palette = "Qualitative.Bold",
    reverse_y = TRUE,
    order = "file",              # "file" or "median"
    outlier_shape = 18,          # 18 = filled diamond
    outlier_size = 2.4,
    labels = list(x = "Subtype", y = "Depth"),
    size = list(width = 10, height = 4)
)

load_packages(c("ggplot2", "readr"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- as.data.frame(read_table_auto(io$input), stringsAsFactors = FALSE)
require_columns(df, config$columns)
df <- data.frame(x = as.character(df[[config$columns$x]]),
    y = as.numeric(df[[config$columns$y]]), stringsAsFactors = FALSE)
df <- df[is.finite(df$y), , drop = FALSE]
if (!nrow(df)) stop("No finite values to plot.", call. = FALSE)

lv <- unique(df$x)
if (identical(config$order, "median")) {
    md <- tapply(df$y, df$x, stats::median)
    lv <- names(sort(md[lv]))
}
cols <- stats::setNames(palette_colors(config$palette, n = length(lv)), lv)

boxes <- list()
meds <- list()
outs <- list()
for (i in seq_along(lv)) {
    v <- df$y[df$x == lv[[i]]]
    n <- length(v)
    k <- if (identical(config$depth, "tukey")) floor(log2(n)) - 3L else
        as.integer(config$depth)
    k <- max(1L, min(k, as.integer(config$max_levels)))
    p_lo <- 0.5^(seq_len(k) + 1)
    q_lo <- stats::quantile(v, p_lo, names = FALSE)
    q_hi <- stats::quantile(v, 1 - p_lo, names = FALSE)
    w <- config$box_width *
        (1 - (seq_len(k) - 1) / max(k, 2) * (1 - config$min_width))
    shade <- grDevices::colorRampPalette(c(cols[[i]], "white"))(k + 2)
    boxes[[i]] <- data.frame(x = lv[[i]], xi = i, level = seq_len(k),
        ymin = q_lo, ymax = q_hi, xmin = i - w / 2, xmax = i + w / 2,
        fill = shade[seq_len(k)])
    meds[[i]] <- data.frame(xi = i, y = stats::median(v),
        x0 = i - config$box_width / 2, x1 = i + config$box_width / 2)
    o <- v[v < min(q_lo) | v > max(q_hi)]
    if (length(o)) outs[[i]] <- data.frame(xi = i, y = o)
}
boxes <- do.call(rbind, boxes)
boxes <- boxes[order(-boxes$level), , drop = FALSE]
meds <- do.call(rbind, meds)
outs <- do.call(rbind, outs)

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p <- ggplot() +
    geom_rect(data = boxes, aes(xmin = xmin, xmax = xmax, ymin = ymin,
        ymax = ymax, fill = fill), colour = "black", linewidth = 0.7) +
    geom_segment(data = meds, aes(x = x0, xend = x1, y = y, yend = y),
        colour = "grey35", linewidth = 0.6) +
    scale_fill_identity() +
    scale_x_continuous(breaks = seq_along(lv), labels = lv,
        expand = expansion(add = 0.5)) +
    labs(x = config$labels$x, y = config$labels$y) +
    theme_classic(base_size = 11) +
    theme(
        axis.text = element_text(colour = "black"),
        axis.ticks.x = element_line(colour = "black"),
        plot.background = element_blank(),
        panel.background = element_blank()
    )
if (!is.null(outs) && nrow(outs)) {
    p <- p + geom_point(data = outs, aes(x = xi, y = y),
        shape = config$outlier_shape, size = config$outlier_size,
        colour = "black")
}
if (isTRUE(config$reverse_y)) {
    p <- p + scale_y_reverse(breaks = scales::pretty_breaks(5))
}

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output, width = config$size$width,
    height = config$size$height)
