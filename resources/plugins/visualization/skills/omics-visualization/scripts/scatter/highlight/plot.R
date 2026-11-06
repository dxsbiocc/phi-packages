#!/usr/bin/env Rscript

# Template-ID: scatter-highlight
#
# Purpose:
#   Small-multiple scatter: one panel per group, the group's points in
#   colour on top of every other point in grey, so each group is read
#   against the whole distribution. Strip colour matches the group.
#
# Inputs:
#   A long table (tsv/csv). Default example (ggplot2::mpg):
#     - displ: x
#     - hwy: y
#     - class: group (one panel each)
#
# Output:
#   A PDF, PNG, or SVG faceted highlight scatter.
#
# Dependencies:
#   ggplot2, readr, ggh4x
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   Edit only CONFIG: columns, labels, palette / group_colours,
#   strip_colour ("group" or "none"), background colour, ncol,
#   free_scales.
#
# Scientific assumptions:
#   Groups are supplied; the script does not cluster. Every panel
#   shares the same axes unless free_scales = TRUE, so positions are
#   comparable across panels. Grey points are all other observations,
#   drawn identically in every panel.

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
        x = "displ",
        y = "hwy",
        group = "class"
    ),
    labels = list(
        x = "displ",
        y = "hwy"
    ),
    palettes = list(group = "Qualitative.Bold"),
    group_colours = NULL,       # named hex per group; NULL = palette
    strip_colour = "group",     # "group" (strip = point colour) or "none"
    strip_alpha = 0.85,
    background_colour = "grey85",
    point_size = 1.8,
    ncol = 3L,
    free_scales = FALSE
)

load_packages(c("ggplot2", "readr", "ggh4x"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- as.data.frame(read_table_auto(io$input), stringsAsFactors = FALSE)
cc <- config$columns
require_columns(df, cc)
df$.x <- as.numeric(df[[cc$x]])
df$.y <- as.numeric(df[[cc$y]])
df$.g <- as.character(df[[cc$group]])
df <- df[is.finite(df$.x) & is.finite(df$.y) & !is.na(df$.g) & nzchar(df$.g), , drop = FALSE]
if (!nrow(df)) {
    stop("No finite x-y rows with a group.", call. = FALSE)
}
g_levels <- unique(df$.g)
if (length(g_levels) < 2L) {
    stop("Need at least two groups.", call. = FALSE)
}
df$.g <- factor(df$.g, levels = g_levels)

gcols <- palette_colors(config$palettes$group, n = length(g_levels))
names(gcols) <- g_levels
if (!is.null(config$group_colours)) {
    hit <- intersect(names(config$group_colours), g_levels)
    gcols[hit] <- config$group_colours[hit]
}

# Background layer: every observation, without the facet variable, so it
# repeats in every panel.
bg <- data.frame(.x = df$.x, .y = df$.y)

text_on <- function(hex) {
    rgb <- grDevices::col2rgb(hex) / 255
    lum <- 0.2126 * rgb[1, ] + 0.7152 * rgb[2, ] + 0.0722 * rgb[3, ]
    ifelse(lum > 0.55, "grey10", "white")
}
strip <- if (identical(config$strip_colour, "group")) {
    ggh4x::strip_themed(
        background_x = lapply(g_levels, function(g) {
            ggplot2::element_rect(
                fill = grDevices::adjustcolor(gcols[[g]], alpha.f = config$strip_alpha),
                colour = "grey20", linewidth = 0.4
            )
        }),
        text_x = lapply(g_levels, function(g) {
            ggplot2::element_text(colour = text_on(gcols[[g]]), size = 11)
        })
    )
} else {
    ggh4x::strip_themed(
        background_x = ggplot2::element_rect(fill = "grey92", colour = "grey20",
            linewidth = 0.4),
        text_x = ggplot2::element_text(colour = "grey10", size = 11)
    )
}

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p <- ggplot2::ggplot(df, ggplot2::aes(x = .x, y = .y)) +
    ggplot2::geom_point(data = bg, colour = config$background_colour,
        size = config$point_size) +
    ggplot2::geom_point(ggplot2::aes(colour = .g), size = config$point_size,
        show.legend = FALSE) +
    ggh4x::facet_wrap2(
        ggplot2::vars(.g),
        ncol = as.integer(config$ncol),
        scales = if (isTRUE(config$free_scales)) "free" else "fixed",
        strip = strip
    ) +
    ggplot2::scale_colour_manual(values = gcols) +
    ggplot2::labs(x = config$labels$x, y = config$labels$y) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        plot.background = ggplot2::element_rect(fill = "white", colour = NA),
        plot.margin = ggplot2::margin(8, 10, 8, 8)
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
n_col <- min(as.integer(config$ncol), length(g_levels))
n_row <- ceiling(length(g_levels) / n_col)
save_ggplot(p, io$output, width = 2.9 * n_col + 0.8, height = 2.1 * n_row + 0.8)
