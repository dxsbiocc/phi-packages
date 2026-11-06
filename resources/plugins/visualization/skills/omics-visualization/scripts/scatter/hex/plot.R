#!/usr/bin/env Rscript

# Template-ID: scatter-hex
#
# Purpose:
#   Draw a hexagonal or rectangular (bin2d) 2D histogram of two
#   numeric variables, with optional lm fit, Pearson label, axis
#   limits, cell borders, and top/right marginal histograms.
#
# Inputs:
#   One row per observation. Default example:
#     - BRCA1, BRCA2
#   Diamonds-style example (example_xy.tsv):
#     - x, y
#
# Output:
#   A PDF, PNG, or SVG 2D bin plot.
#
# Dependencies:
#   ggplot2, readr, ggpubr, hexbin, patchwork
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   Square bins (bin_shape), bins / binwidth, fit on/off, cell
#   border, xlim/ylim and marginal histograms are CONFIG, not new
#   ids. 2D kernel densities stay scatter-density2d.
#
# Scientific assumptions:
#   Each row is one paired observation. Bin counts are a display
#   binning of the raw points. The optional lm / Pearson label is a
#   display fit, not a causal claim.

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
    columns = list(
        x = "BRCA1",
        y = "BRCA2"
    ),
    bin_shape = "hex",           # "hex" or "square" (geom_bin2d)
    bins = 30,
    binwidth = NULL,             # c(wx, wy) overrides bins
    show_fit = TRUE,             # lm line + Pearson label
    marginal = "none",           # "none" or "histogram"
    marginal_bins = 35,
    xlim = NULL,                 # e.g. c(4, 10)
    ylim = NULL,
    cell_colour = "grey70",
    cell_linewidth = 0.2,
    fill_palette = "Gradient.YlGnBu",
    labels = list(
        title = "Correlation of BRCA1 and BRCA2",
        x = "BRCA1",
        y = "BRCA2"
    ),
    size = list(width = 7.5, height = 6.5)
)

load_packages(c("ggplot2", "readr", "ggpubr", "hexbin", "patchwork"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- as.data.frame(read_table_auto(io$input), stringsAsFactors = FALSE)
require_columns(df, config$columns)
x_col <- config$columns$x
y_col <- config$columns$y
df[[x_col]] <- as.numeric(df[[x_col]])
df[[y_col]] <- as.numeric(df[[y_col]])
ok <- is.finite(df[[x_col]]) & is.finite(df[[y_col]])
df <- df[ok, , drop = FALSE]
if (!is.null(config$xlim)) {
    df <- df[df[[x_col]] >= config$xlim[[1]] & df[[x_col]] <= config$xlim[[2]], , drop = FALSE]
}
if (!is.null(config$ylim)) {
    df <- df[df[[y_col]] >= config$ylim[[1]] & df[[y_col]] <= config$ylim[[2]], , drop = FALSE]
}
if (nrow(df) < 5L) stop("Too few points after filtering.", call. = FALSE)

fill_stops <- palette_colors(config$fill_palette)

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
bin_args <- list(
    colour = config$cell_colour,
    linewidth = config$cell_linewidth
)
if (is.null(config$binwidth)) {
    bin_args$bins <- config$bins
} else {
    bin_args$binwidth <- config$binwidth
}
bin_geom <- if (identical(config$bin_shape, "square")) {
    ggplot2::geom_bin2d
} else {
    ggplot2::geom_hex
}

p <- ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    do.call(bin_geom, bin_args) +
    scale_fill_gradientn(colours = fill_stops, name = "count")
if (isTRUE(config$show_fit)) {
    p <- p +
        geom_smooth(
            method = "lm", formula = y ~ x,
            colour = palette_colors("Qualitative.Bold")[[8]],
            fill = "grey80", alpha = 0.25, linewidth = 0.7
        ) +
        ggpubr::stat_cor(method = "pearson", size = 4.5)
}
if (!is.null(config$xlim) || !is.null(config$ylim)) {
    p <- p + coord_cartesian(xlim = config$xlim, ylim = config$ylim)
}
p <- p +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y
    ) +
    theme_bw(base_size = 13) +
    theme(
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(
            colour = "grey88", linetype = "dashed", linewidth = 0.3
        ),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
        axis.text = element_text(colour = "black"),
        plot.background = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA),
        legend.background = element_rect(fill = "white", colour = NA)
    )
if (is.null(config$labels$title) || !nzchar(as.character(config$labels$title))) {
    p <- p + theme(plot.title = element_blank())
}

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
if (identical(config$marginal, "histogram")) {
    hist_fill <- fill_stops[[max(1L, length(fill_stops) - 2L)]]
    void <- theme_void() + theme(plot.background = element_blank())
    top <- ggplot(df, aes(x = .data[[x_col]])) +
        geom_histogram(
            bins = config$marginal_bins, fill = hist_fill,
            colour = "grey20", linewidth = 0.2
        ) + void
    right <- ggplot(df, aes(y = .data[[y_col]])) +
        geom_histogram(
            bins = config$marginal_bins, fill = hist_fill,
            colour = "grey20", linewidth = 0.2
        ) + void
    p <- top + patchwork::plot_spacer() + (p + labs(title = NULL)) + right +
        patchwork::plot_layout(
            ncol = 2, widths = c(4, 1), heights = c(1, 4),
            guides = "collect"
        ) +
        patchwork::plot_annotation(
            title = config$labels$title,
            theme = theme(
                plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
            )
        )
    save_ggplot(
        p, io$output,
        width = config$size$width + 0.5,
        height = config$size$height + 0.8
    )
} else {
    save_ggplot(
        p, io$output,
        width = config$size$width,
        height = config$size$height
    )
}
