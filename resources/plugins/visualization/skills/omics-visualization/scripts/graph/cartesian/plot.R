#!/usr/bin/env Rscript

# Template-ID: graph-cartesian
#
# Purpose:
#   Draw a path graph of values on a cartesian category axis.
#
# Inputs:
#   Two tables (tsv/csv):
#     1. nodes: name, y
#     2. links: source, target (0-based node indices)
#
# Output:
#   A PDF, PNG, or SVG cartesian graph.
#
# Dependencies:
#   ggplot2, readr, ggprism
#
# Example:
#   Rscript plot.R nodes.tsv links.tsv output.pdf
#
# Agent adaptation:
#   For a new dataset, edit only CONFIG (column names and labels).
#   Pass the node table then the edge table on the CLI.
#   Edit PLOT only when the chart geometry or styling must change.
#
# Scientific assumptions:
#   Nodes are ordered as in the file; links use 0-based row indices.
#   The path encodes a supplied edge list, not a computed trajectory.
#   Y values are already summarized observations.

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

io <- parse_io_args(n_input = 2L, input_names = c("nodes", "links"))

# -----------------------------------------------------------------------------
# CONFIG  (edit this block for a new dataset)
# -----------------------------------------------------------------------------
config <- list(
    columns = list(
        x = "name",
        y = "y"
    ),
    labels = list(
        title = "Cartesian Graph",
        x = "Day",
        y = "Value"
    )
)

load_packages(c("ggplot2", "readr", "ggprism"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
nodes <- read_table_auto(io$input[[1]])
require_columns(nodes, config$columns)
x_col <- config$columns$x
y_col <- config$columns$y
nodes[[x_col]] <- factor(nodes[[x_col]], levels = unique(nodes[[x_col]]))
nodes[[y_col]] <- as.numeric(nodes[[y_col]])
nodes$x_num <- as.numeric(nodes[[x_col]])

links <- read_table_auto(io$input[[2]])
require_columns(links, list(source = "source", target = "target"))
from_i <- as.integer(links$source) + 1L
to_i <- as.integer(links$target) + 1L
edge_df <- data.frame(
    x = nodes$x_num[from_i],
    y = nodes[[y_col]][from_i],
    xend = nodes$x_num[to_i],
    yend = nodes[[y_col]][to_i]
)

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p <- ggplot() +
    geom_segment(
        data = edge_df,
        aes(x = x, y = y, xend = xend, yend = yend),
        colour = "#2f4554",
        linewidth = 1
    ) +
    geom_point(
        data = nodes,
        aes(x = x_num, y = .data[[y_col]]),
        size = 3,
        colour = "#1cc7d0"
    ) +
    scale_x_continuous(breaks = nodes$x_num, labels = as.character(nodes[[x_col]])) +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y
    ) +
    theme_prism() +
    theme(plot.background = element_blank())

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output)
