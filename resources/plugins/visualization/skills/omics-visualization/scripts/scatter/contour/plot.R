#!/usr/bin/env Rscript

# Template-ID: scatter-contour
#
# Purpose:
#   Draw a scatter plot of two numeric axes with a convex-hull outline
#   around each cluster.
#
# Inputs:
#   A table with one row per point. Default example:
#     - x: first numeric axis
#     - y: second numeric axis
#     - cluster: group
#
# Output:
#   A PDF, PNG, or SVG contour scatter plot.
#
# Dependencies:
#   ggplot2, readr
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (column names and labels).
#   Edit DATA PREPARATION to change cluster order or how hulls are
#   computed.
#   Edit PLOT only when the chart geometry or styling must change.
#
# Scientific assumptions:
#   Each row is one observation with two coordinates and a cluster label.
#   Convex hulls are a display of the cluster extent, not a density
#   estimate or statistical contour.

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
    labels = list(
        title = "Scatter Plot with Contours",
        x = "pc1",
        y = "pc2"
    )
)

load_packages(c("ggplot2", "readr"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
require_columns(df, config$columns)
x_col <- config$columns$x
y_col <- config$columns$y
fill_col <- config$columns$fill
df[[fill_col]] <- factor(df[[fill_col]], levels = unique(df[[fill_col]]))
df[[x_col]] <- as.numeric(df[[x_col]])
df[[y_col]] <- as.numeric(df[[y_col]])

hull_data <- do.call(rbind, lapply(
    split(df, df[[fill_col]]),
    function(sub) {
        sub[grDevices::chull(sub[[x_col]], sub[[y_col]]), , drop = FALSE]
    }
))

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
fill_values <- c("#1CC7D0", "#2DDE98", "#FFC168")
colour_values <- c("#1cc7d0", "#2dde98", "#ffc168")

p <- ggplot(df, aes(
    x = .data[[x_col]],
    y = .data[[y_col]],
    fill = .data[[fill_col]]
)) +
    geom_vline(xintercept = 0, linewidth = 0.2) +
    geom_hline(yintercept = 0, linewidth = 0.2) +
    geom_point(shape = 21, size = 4, alpha = 0.8) +
    geom_polygon(
        data = hull_data,
        aes(
            x = .data[[x_col]],
            y = .data[[y_col]],
            fill = .data[[fill_col]],
            colour = .data[[fill_col]]
        ),
        alpha = 0.2
    ) +
    scale_fill_manual(values = fill_values) +
    scale_colour_manual(values = colour_values) +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y
    ) +
    theme_minimal()

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output)
