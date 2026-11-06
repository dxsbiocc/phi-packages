#!/usr/bin/env Rscript

# Template-ID: line-ribbon-diff
#
# Purpose:
#   Draw two numeric series on a shared ordered x axis, fill the
#   band between them with one colour when series A is above B and
#   another when B is above A, and overlay both lines.
#
# Inputs:
#   Wide table, one row per x. Default example (economics pce /
#   unemploy, value01-scaled):
#     - date, pce, unemploy
#
# Output:
#   A PDF, PNG, or SVG difference-ribbon line chart.
#
# Dependencies:
#   ggplot2, readr
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   Edit CONFIG for column names, palette refs, and labels.
#   Dual-axis hydrographs stay line-double. Plain multi-series
#   overlays stay line-group. Do not rescale series beyond the
#   supplied values; value01-style scaling is upstream.
#
# Scientific assumptions:
#   Each row is one shared x with two already measured series.
#   Ribbon colours mark which series is higher at that x; they are
#   a display of the supplied values, not a test.

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
        x = "date",
        a = "pce",
        b = "unemploy"
    ),
    labels = list(
        title = NULL,
        x = NULL,
        y = NULL,
        a = "pce",
        b = "unemploy",
        a_above = "pce > unemploy",
        b_above = "unemploy > pce"
    ),
    # Catalog colours matching the classic Set3 economics ribbon example.
    palettes = list(
        ribbon_a_above = "Qualitative.Set3 / 12:6",  # #fdb462
        ribbon_b_above = "Qualitative.Set3 / 12:1",  # #8dd3c7
        line_a = "Qualitative.Set3 / 12:4",          # #fb8072
        line_b = "Qualitative.Set3 / 12:5"           # #80b1d3
    ),
    ribbon_alpha = 0.7,
    line_width = 0.75,
    x_is_date = TRUE,
    show_legend = TRUE,
    base_size = 13,
    size = list(width = 9, height = 4.8)
)

load_packages(c("ggplot2", "readr"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
pal_ref <- function(ref) {
    parts <- strsplit(ref, ":", fixed = TRUE)[[1]]
    cols <- palette_colors(parts[[1]])
    if (length(parts) > 1L) cols[[as.integer(parts[[2]])]] else cols
}

df <- as.data.frame(read_table_auto(io$input), stringsAsFactors = FALSE)
require_columns(df, config$columns)
x_col <- config$columns$x
a_col <- config$columns$a
b_col <- config$columns$b

if (isTRUE(config$x_is_date)) {
    df[[x_col]] <- as.Date(df[[x_col]])
} else {
    df[[x_col]] <- as.numeric(df[[x_col]])
}
df[[a_col]] <- as.numeric(df[[a_col]])
df[[b_col]] <- as.numeric(df[[b_col]])
ok <- !is.na(df[[x_col]]) & is.finite(df[[a_col]]) & is.finite(df[[b_col]])
if (any(!ok)) {
    message("Dropped ", sum(!ok), " incomplete row(s).")
    df <- df[ok, , drop = FALSE]
}
if (nrow(df) < 2L) stop("Need at least 2 rows.", call. = FALSE)
df <- df[order(df[[x_col]]), , drop = FALSE]

df$low_min <- pmin(df[[a_col]], df[[b_col]])
df$low_max <- pmax(df[[a_col]], df[[b_col]])
df$high_min <- df$low_min
df$high_max <- df$low_max
# low_* : when A <= B (B above); high_* : when A > B (A above)
a_above <- df[[a_col]] > df[[b_col]]
df$low_min[a_above] <- NA
df$low_max[a_above] <- NA
df$high_min[!a_above] <- NA
df$high_max[!a_above] <- NA

col_ribbon_a <- pal_ref(config$palettes$ribbon_a_above)
col_ribbon_b <- pal_ref(config$palettes$ribbon_b_above)
col_line_a <- pal_ref(config$palettes$line_a)
col_line_b <- pal_ref(config$palettes$line_b)

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p <- ggplot(df, aes(x = .data[[x_col]])) +
    geom_ribbon(
        aes(ymin = low_min, ymax = low_max, fill = "b_above"),
        alpha = config$ribbon_alpha, na.rm = TRUE
    ) +
    geom_ribbon(
        aes(ymin = high_min, ymax = high_max, fill = "a_above"),
        alpha = config$ribbon_alpha, na.rm = TRUE
    ) +
    geom_line(
        aes(y = .data[[a_col]], colour = "a"),
        linewidth = config$line_width
    ) +
    geom_line(
        aes(y = .data[[b_col]], colour = "b"),
        linewidth = config$line_width
    ) +
    scale_fill_manual(
        name = NULL,
        values = c(a_above = col_ribbon_a, b_above = col_ribbon_b),
        labels = c(
            a_above = config$labels$a_above,
            b_above = config$labels$b_above
        ),
        breaks = c("a_above", "b_above"),
        guide = if (isTRUE(config$show_legend)) "legend" else "none"
    ) +
    scale_colour_manual(
        name = NULL,
        values = c(a = col_line_a, b = col_line_b),
        labels = c(a = config$labels$a, b = config$labels$b),
        breaks = c("a", "b"),
        guide = if (isTRUE(config$show_legend)) "legend" else "none"
    ) +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y
    ) +
    theme_classic(base_size = config$base_size) +
    theme(
        axis.text = element_text(colour = "black"),
        legend.position = "bottom",
        legend.box = "horizontal",
        plot.background = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA)
    )

if (is.null(config$labels$title) || !nzchar(as.character(config$labels$title))) {
    p <- p + theme(plot.title = element_blank())
}
if (is.null(config$labels$x) || !nzchar(as.character(config$labels$x))) {
    p <- p + theme(axis.title.x = element_blank())
}
if (is.null(config$labels$y) || !nzchar(as.character(config$labels$y))) {
    p <- p + theme(axis.title.y = element_blank())
}

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(
    p, io$output,
    width = config$size$width,
    height = config$size$height
)
