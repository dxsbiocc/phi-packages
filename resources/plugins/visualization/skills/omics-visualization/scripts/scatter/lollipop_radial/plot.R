#!/usr/bin/env Rscript

# Template-ID: scatter-lollipop-radial
#
# Purpose:
#   Draw a radial lollipop chart of a numeric value around a categorical
#   axis, using a fan-shaped polar coordinate system.
#
# Inputs:
#   A table with one row per category. Default example:
#     - age_group: category label
#     - mean_prev: numeric value
#     - group: colour group
#
# Output:
#   A PDF, PNG, or SVG radial lollipop plot.
#
# Dependencies:
#   ggplot2, readr
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (column names and labels).
#   Edit DATA PREPARATION to change category order.
#   Edit PLOT only when the chart geometry or styling must change.
#
# Scientific assumptions:
#   Each row is one category with a comparable numeric value.
#   Polar layout is a display choice, not a circular data transform.
#   Values are treated as already summarized.

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
        x = "age_group",
        y = "mean_prev",
        colour = "group"
    ),
    labels = list(
        title = "Radial Lollipop Plot",
        x = "",
        y = ""
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
colour_col <- config$columns$colour

df[[x_col]] <- factor(df[[x_col]], levels = unique(df[[x_col]]))
df[[colour_col]] <- factor(df[[colour_col]], levels = unique(df[[colour_col]]))
df[[y_col]] <- as.numeric(df[[y_col]])

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
colour_values <- c("#1aafd0", "#2dde98", "#ffc168")

p <- ggplot(df, aes(
    x = .data[[config$columns$x]],
    y = .data[[config$columns$y]],
    colour = .data[[config$columns$colour]]
)) +
    geom_segment(aes(
        x = .data[[config$columns$x]],
        xend = .data[[config$columns$x]],
        y = 0,
        yend = .data[[config$columns$y]]
    )) +
    geom_point(size = 5, shape = 21, fill = "#ffffff") +
    geom_point(size = 3) +
    geom_hline(yintercept = 0, colour = "#81947a") +
    coord_radial(start = -1.256637, end = 1.256637, inner.radius = 0.4) +
    scale_colour_manual(values = colour_values) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    guides(colour = "none") +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y
    ) +
    theme(
        panel.background = element_rect(fill = "#f3ead833"),
        panel.grid.major = element_line(colour = "grey90"),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 14, hjust = 0.5),
        plot.background = element_blank()
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output)
