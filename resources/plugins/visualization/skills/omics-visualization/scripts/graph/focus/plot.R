#!/usr/bin/env Rscript

# Template-ID: graph-focus
#
# Purpose:
#   Place one named node at the origin and arrange the rest by graph
#   distance on concentric rings.
#
# Inputs:
#   An edge table with one row per link. Default example:
#     - from: source node
#     - to: target node
#     - corr: signed correlation or similar strength in [-1, 1]
#   Set config$focus to the node that should sit at the center.
#
# Output:
#   A PDF, PNG, or SVG focus network plot.
#
# Dependencies:
#   ggplot2, readr, ggraph, tidygraph, ggprism, ggforce
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (column names, focus name, labels).
#   Edit DATA PREPARATION to change directed vs undirected layout.
#   Edit PLOT only when the chart geometry or styling must change.
#
# Scientific assumptions:
#   Focus identity is supplied, not chosen by centrality.
#   Ring radius is graph distance from that node, a layout coordinate,
#   not a statistical neighborhood test.
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
    focus = "TP53",
    labels = list(
        title = "Focus network",
        x = "",
        y = "",
        edge = "Sign",
        width = "|corr|"
    )
)

load_packages(c("ggplot2", "readr", "ggraph", "tidygraph", "ggprism", "ggforce"))

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
if (!config$focus %in% node_order) {
    stop(
        paste("Focus node is not in the edge table:", config$focus),
        call. = FALSE
    )
}
nodes <- data.frame(
    name = node_order,
    is_focus = node_order == config$focus,
    stringsAsFactors = FALSE
)
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
    directed = FALSE
)
focus_idx <- match(config$focus, nodes$name)
layout <- ggraph::create_layout(
    graph,
    layout = "focus",
    focus = focus_idx
)
ring_n <- max(as.integer(layout$distance), na.rm = TRUE)
if (!is.finite(ring_n) || ring_n < 1L) {
    ring_n <- 1L
}
rings <- data.frame(r = seq_len(ring_n), stringsAsFactors = FALSE)

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
title <- config$labels$title
if (!nzchar(title)) {
    title <- paste("Focus:", config$focus)
}

p <- ggraph::ggraph(layout) +
    ggforce::geom_circle(
        data = rings,
        ggplot2::aes(x0 = 0, y0 = 0, r = r),
        inherit.aes = FALSE,
        colour = "#d8dce3",
        linewidth = 0.5
    ) +
    ggraph::geom_edge_link(
        ggplot2::aes(colour = corr_label, width = corr_abs),
        alpha = 0.9
    ) +
    ggraph::geom_node_point(
        ggplot2::aes(size = is_focus, colour = is_focus)
    ) +
    ggraph::geom_node_text(
        ggplot2::aes(
            label = name,
            colour = is_focus,
            fontface = ifelse(is_focus, "bold", "plain")
        ),
        size = 4,
        repel = FALSE,
        vjust = -1
    ) +
    ggplot2::scale_size_manual(
        values = c("FALSE" = 3.6, "TRUE" = 6.5),
        guide = "none"
    ) +
    ggplot2::scale_colour_manual(
        values = c("FALSE" = "#1aafd0", "TRUE" = "#6a67ce"),
        guide = "none"
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
        range = c(0.45, 2.8),
        limits = c(0, 1)
    ) +
    ggplot2::scale_x_continuous(expand = expansion(mult = 0.12)) +
    ggplot2::scale_y_continuous(expand = expansion(mult = 0.12)) +
    ggplot2::coord_fixed(clip = "off") +
    ggplot2::labs(
        title = title,
        subtitle = paste("Center:", config$focus),
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
        plot.subtitle = element_text(size = 11),
        plot.background = element_blank(),
        plot.margin = margin(12, 12, 12, 12)
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output, width = 8, height = 8)
