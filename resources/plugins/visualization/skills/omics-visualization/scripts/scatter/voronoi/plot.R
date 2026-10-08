#!/usr/bin/env Rscript

# Template-ID: scatter-voronoi
#
# Purpose:
#   Draw a Voronoi tessellation of supplied x-y points, filled by a
#   supplied group label, with the points overlaid.
#
# Inputs:
#   One row per observation. Default example:
#     - x, y: coordinates (for example a supplied embedding)
#     - cluster: group label
#
# Output:
#   A PDF, PNG, or SVG Voronoi scatter plot.
#
# Dependencies:
#   ggplot2, readr, ggprism, ggforce
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (columns, max.radius, expand,
#   point_size).
#   Edit DATA PREPARATION to change group order.
#   Edit PLOT only when the chart geometry or styling must change.
#
# Scientific assumptions:
#   x, y, and group are supplied. This script does not run PCA, UMAP,
#   or clustering. Tile area is a nearest-site partition of the plane,
#   not a density or a statistical region.

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
        x = "x",
        y = "y",
        fill = "cluster"
    ),
    palette = "Brand.Algolia",
    # NULL = tiles fill the panel. A positive number clips each cell.
    max_radius = NULL,
    expand_mm = -0.45,
    radius_mm = 1.2,
    point_size = 3.2,
    labels = list(
        title = "",
        x = "Dim 1",
        y = "Dim 2",
        fill = NULL
    )
)

load_packages(c("ggplot2", "readr", "ggprism", "ggforce"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
require_columns(df, config$columns)

x_col <- config$columns$x
y_col <- config$columns$y
fill_col <- config$columns$fill

df[[x_col]] <- as.numeric(df[[x_col]])
df[[y_col]] <- as.numeric(df[[y_col]])
ok <- is.finite(df[[x_col]]) & is.finite(df[[y_col]])
if (any(!ok)) {
    message("Dropped ", sum(!ok), " row(s) with non-finite x or y.")
    df <- df[ok, , drop = FALSE]
}
if (!nrow(df)) {
    stop("No rows left to plot.", call. = FALSE)
}
df[[fill_col]] <- factor(df[[fill_col]], levels = unique(df[[fill_col]]))

fill_n <- nlevels(df[[fill_col]])
fill_cols <- palette_colors(config$palette, n = fill_n)
border_col <- "#ffffff"
point_col <- palette_colors("Qualitative.Bold")[[11]]
max_radius <- config$max_radius
if (!is.null(max_radius) && (!is.finite(max_radius) || max_radius <= 0)) {
    max_radius <- NULL
}
fill_name <- config$labels$fill
if (is.null(fill_name)) {
    fill_name <- fill_col
}

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p <- ggplot(df, aes(
    x = .data[[x_col]],
    y = .data[[y_col]],
    group = -1L
)) +
    geom_voronoi_tile(
        aes(fill = .data[[fill_col]]),
        colour = border_col,
        linewidth = 0.25,
        max.radius = max_radius,
        expand = unit(config$expand_mm, "mm"),
        radius = unit(config$radius_mm, "mm"),
        show.legend = TRUE
    ) +
    geom_point(
        shape = 21,
        size = config$point_size,
        fill = border_col,
        colour = point_col,
        stroke = 0.55,
        show.legend = FALSE
    ) +
    scale_fill_manual(values = fill_cols, name = fill_name) +
    coord_fixed(clip = "off") +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y
    ) +
    theme_prism() +
    theme(
        panel.grid = element_blank(),
        axis.line = element_line(colour = "#888888", linewidth = 0.6),
        axis.ticks = element_line(colour = "#888888"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        legend.background = element_blank()
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output, width = 7.5, height = 6.2)
