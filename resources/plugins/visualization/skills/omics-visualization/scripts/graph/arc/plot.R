#!/usr/bin/env Rscript

# Template-ID: graph-arc
#
# Purpose:
#   Place nodes on a line and draw signed edges as arcs above the axis.
#
# Inputs:
#   An edge table with one row per link. Default example:
#     - from: source node
#     - to: target node
#     - corr: signed correlation or similar strength in [-1, 1]
#
# Output:
#   A PDF, PNG, or SVG arc diagram.
#
# Dependencies:
#   ggplot2, readr, ggraph, tidygraph, ggprism
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (column names and labels).
#   Edit DATA PREPARATION to change node order along the line.
#   Edit PLOT only when the chart geometry or styling must change.
#
# Scientific assumptions:
#   Linear order is a display choice, not a clustering or trajectory result.
#   corr is a supplied signed strength. Colour encodes sign; edge width
#   encodes absolute magnitude. This script does not compute correlation.

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
        from = "from",
        to = "to",
        corr = "corr"
    ),
    labels = list(
        title = "Arc diagram",
        x = "",
        y = "",
        edge = "Sign",
        width = "|corr|"
    )
)

load_packages(c("ggplot2", "readr", "ggraph", "tidygraph", "ggprism"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
edges <- read_table_auto(io$input)
require_columns(edges, config$columns)
from_col <- config$columns$from
to_col <- config$columns$to
corr_col <- config$columns$corr
edges[[from_col]] <- as.character(edges[[from_col]])
edges[[to_col]] <- as.character(edges[[to_col]])
edges[[corr_col]] <- as.numeric(edges[[corr_col]])
if (any(!is.finite(edges[[corr_col]]))) {
    stop("corr must be numeric and finite.", call. = FALSE)
}
edges$corr_abs <- abs(edges[[corr_col]])
edges$corr_label <- factor(
    ifelse(
        edges[[corr_col]] > 0,
        "Positive",
        ifelse(edges[[corr_col]] < 0, "Negative", "Neutral")
    ),
    levels = c("Negative", "Neutral", "Positive")
)
node_order <- unique(c(edges[[from_col]], edges[[to_col]]))
nodes <- data.frame(name = node_order, stringsAsFactors = FALSE)
graph <- tidygraph::tbl_graph(
    nodes = nodes,
    edges = data.frame(
        from = edges[[from_col]],
        to = edges[[to_col]],
        corr_label = edges$corr_label,
        corr_abs = edges$corr_abs,
        stringsAsFactors = FALSE
    ),
    node_key = "name",
    directed = TRUE
)

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p <- ggraph::ggraph(graph, layout = "linear") +
    ggraph::geom_edge_arc(
        ggplot2::aes(colour = corr_label, width = corr_abs),
        strength = 0.8,
        alpha = 0.88,
        fold = TRUE
    ) +
    ggraph::geom_node_point(size = 3.2, colour = "#333333") +
    ggraph::geom_node_text(
        ggplot2::aes(label = name),
        angle = 90,
        hjust = 1,
        nudge_y = -0.08,
        size = 4,
        colour = "#333333"
    ) +
    ggraph::scale_edge_colour_manual(
        name = config$labels$edge,
        values = c(
            Negative = "#fc636b",
            Neutral = "#a8adb7",
            Positive = "#3be8b0"
        )
    ) +
    ggraph::scale_edge_width(
        name = config$labels$width,
        range = c(0.4, 2.4),
        limits = c(0, 1)
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y
    ) +
    theme_prism(base_size = 12) +
    theme(
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 15),
        plot.background = element_blank(),
        plot.margin = margin(12, 16, 80, 16)
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output, width = 14, height = 6)
