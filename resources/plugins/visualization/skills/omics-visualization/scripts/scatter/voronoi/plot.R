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
#   ggplot2, readr, ggprism, ggforce, deldir
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (columns, max.radius, bound,
#   normalize, show_edges, overlay, expand, point_size). A radius
#   clip, a hull / supplied outline bound and a Delaunay overlay are
#   variants (preview_radius / preview_bound / preview_delaunay), not
#   new ids.
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
    # NULL = rectangle around the points; "hull" = convex hull of the
    # points; a filename = supplied polygon sidecar (x, y) beside the
    # input, e.g. a tissue outline.
    bound = NULL,
    normalize = FALSE,           # tessellate on unit-scaled x / y
    show_edges = FALSE,          # overlay Voronoi cell edges
    # "none", "delaunay" (triangulation over the tiles) or
    # "delaunay_only" (triangulation and group-filled points, no tiles)
    overlay = "none",
    delaunay_max_edge = NULL,    # drop triangulation edges longer than this
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

load_packages(c("ggplot2", "readr", "ggprism", "ggforce", "deldir"))

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

bound <- NULL
if (identical(config$bound, "hull")) {
    h <- grDevices::chull(df[[x_col]], df[[y_col]])
    bound <- cbind(df[[x_col]][h], df[[y_col]][h])
} else if (is.character(config$bound) && nzchar(config$bound)) {
    b_path <- file.path(dirname(io$input), config$bound)
    if (!file.exists(b_path)) {
        stop(paste("Missing bound sidecar:", b_path), call. = FALSE)
    }
    bd <- read_table_auto(b_path)
    require_columns(bd, list(x = "x", y = "y"))
    bound <- cbind(as.numeric(bd$x), as.numeric(bd$y))
    if (nrow(bound) < 3L || any(!is.finite(bound))) {
        stop("Bound polygon needs >= 3 finite vertices.", call. = FALSE)
    }
}
overlay <- config$overlay
draw_tiles <- !identical(overlay, "delaunay_only")

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p <- ggplot(df, aes(
    x = .data[[x_col]],
    y = .data[[y_col]],
    group = -1L
))
if (draw_tiles) {
    p <- p + geom_voronoi_tile(
        aes(fill = .data[[fill_col]]),
        colour = border_col,
        linewidth = 0.25,
        max.radius = max_radius,
        bound = bound,
        normalize = isTRUE(config$normalize),
        expand = unit(config$expand_mm, "mm"),
        radius = unit(config$radius_mm, "mm"),
        show.legend = TRUE
    )
}
if (isTRUE(config$show_edges)) {
    p <- p + geom_voronoi_segment(bound = bound,
        normalize = isTRUE(config$normalize), colour = point_col,
        linewidth = 0.3)
}
if (overlay %in% c("delaunay", "delaunay_only")) {
    tri <- deldir::deldir(df[[x_col]], df[[y_col]])$delsgs
    tri$len <- sqrt((tri$x2 - tri$x1)^2 + (tri$y2 - tri$y1)^2)
    if (!is.null(config$delaunay_max_edge)) {
        tri <- tri[tri$len <= config$delaunay_max_edge, , drop = FALSE]
    }
    p <- p + geom_segment(data = tri, aes(x = x1, y = y1, xend = x2,
        yend = y2), inherit.aes = FALSE, colour = point_col,
        linewidth = 0.4, alpha = if (draw_tiles) 0.75 else 0.9)
}
if (!is.null(bound) && draw_tiles) {
    bd_df <- data.frame(x = bound[, 1], y = bound[, 2])
    p <- p + geom_polygon(data = bd_df, aes(x = x, y = y), inherit.aes = FALSE,
        fill = NA, colour = point_col, linewidth = 0.5, linetype = "dashed")
}
p <- p +
    geom_point(
        if (draw_tiles) NULL else aes(fill = .data[[fill_col]]),
        shape = 21,
        size = config$point_size,
        fill = if (draw_tiles) border_col else NULL,
        colour = point_col,
        stroke = 0.55,
        show.legend = !draw_tiles
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
