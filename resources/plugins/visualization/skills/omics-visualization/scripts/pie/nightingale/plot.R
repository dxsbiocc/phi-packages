#!/usr/bin/env Rscript

# Template-ID: pie-nightingale
#
# Purpose:
#   Draw a Nightingale rose: equal-angle polar bars, radius encodes
#   value, with rounded caps and an inner hole.
#
# Inputs:
#   A two-column table (tsv/csv). Default example:
#     - roseName: category label
#     - value: numeric value
#
# Output:
#   A PDF, PNG, or SVG Nightingale rose chart.
#
# Dependencies:
#   ggplot2, readr, ggprism, gground, ggrepel, grid
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (columns, inner hole, corner
#   radius, value labels). Rounded caps vs a sharp rose, and the
#   inner hole, stay on this id. That is not pie-doughnut-round
#   (angle encodes share of a whole). Edit DATA PREPARATION to
#   change category order. Edit PLOT only when the glyph must change.
#
# Scientific assumptions:
#   Each row is one category with a comparable numeric value.
#   Equal angular width with radius encoding value is a display
#   choice. The inner hole and corner radius are display offsets,
#   not missing data. Values are already summarized.

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
        x = "roseName",
        y = "value"
    ),
    bar_width = 0.95,
    corner_pt = 5,
    inner_frac = 0.28,
    label_frac = 0.28,
    show_value_labels = TRUE,
    labels = list(
        title = "",
        x = "",
        y = "",
        fill = NULL
    )
)

load_packages(c("ggplot2", "readr", "ggprism", "gground", "ggrepel"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
require_columns(df, config$columns)
x_col <- config$columns$x
y_col <- config$columns$y
df[[x_col]] <- factor(df[[x_col]], levels = unique(df[[x_col]]))
df[[y_col]] <- as.numeric(df[[y_col]])
ok <- is.finite(df[[y_col]]) & df[[y_col]] >= 0
if (any(!ok)) {
    message("Dropped ", sum(!ok), " row(s) with invalid values.")
    df <- df[ok, , drop = FALSE]
}
if (!nrow(df)) {
    stop("No rows left to plot.", call. = FALSE)
}

y_max <- max(df[[y_col]], na.rm = TRUE)
if (!is.finite(y_max) || y_max <= 0) {
    stop("Need a positive numeric value for the rose radius.", call. = FALSE)
}
inner_frac <- as.numeric(config$inner_frac)
if (!is.finite(inner_frac) || inner_frac < 0) {
    inner_frac <- 0.28
}
label_frac <- as.numeric(config$label_frac)
if (!is.finite(label_frac) || label_frac < 0) {
    label_frac <- 0.28
}
inner <- y_max * inner_frac
y_top <- y_max + y_max * label_frac
bar_width <- as.numeric(config$bar_width)
if (!is.finite(bar_width) || bar_width <= 0 || bar_width > 1) {
    bar_width <- 0.95
}
corner_pt <- as.numeric(config$corner_pt)
if (!is.finite(corner_pt) || corner_pt < 0) {
    corner_pt <- 5
}

# Reference rose fills from the source figure (圆角扇形图).
rose_seq <- c(
    "#5070dd", "#b6d634", "#ff994d", "#0ca8df",
    "#ffd10a", "#fb628b", "#785db0", "#3fbe95"
)
n_fill <- nlevels(df[[x_col]])
if (n_fill <= length(rose_seq)) {
    fill_vals <- rose_seq[seq_len(n_fill)]
} else {
    fill_vals <- palette_colors("Brand.Algolia", n = n_fill)
    fill_vals[seq_along(rose_seq)] <- rose_seq
}
ink <- palette_colors("Brand.Emma")[[1]]
seg <- palette_colors("Qualitative.Safe")[[12]]

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p <- ggplot(df, aes(
    x = .data[[x_col]],
    y = .data[[y_col]],
    fill = .data[[x_col]]
)) +
    geom_round_col(
        width = bar_width,
        radius = grid::unit(corner_pt, "pt"),
        colour = NA
    ) +
    scale_fill_manual(name = config$labels$fill, values = fill_vals) +
    coord_polar(clip = "off") +
    scale_y_continuous(
        limits = c(-inner, y_top),
        expand = c(0, 0)
    ) +
    labs(
        title = config$labels$title,
        x = NULL,
        y = NULL
    ) +
    guides(fill = guide_legend(nrow = 1, title = NULL)) +
    theme_prism() +
    theme(
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        panel.grid = element_blank(),
        plot.background = element_blank(),
        legend.background = element_blank(),
        legend.position = "bottom",
        legend.key = element_blank()
    )

if (isTRUE(config$show_value_labels)) {
    p <- p + geom_text_repel(
        aes(label = .data[[y_col]]),
        size = 3.6,
        nudge_y = y_max * 0.22,
        segment.curvature = -0.1,
        segment.ncp = 3,
        segment.angle = 20,
        segment.size = 0.4,
        segment.colour = seg,
        colour = ink,
        min.segment.length = 0,
        show.legend = FALSE
    )
}

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output, width = 8, height = 7.2)
