#!/usr/bin/env Rscript

# Template-ID: graph-stress
#
# Purpose:
#   Place nodes with a stress-majorization layout computed from the edge
#   table. Coordinates are not supplied. Do not add separate kk / fr / lgl
#   templates; change layout = "stress" in PLOT only if a sibling algorithm
#   is required.
#
# Inputs:
#   An edge table with one row per link. Default example:
#     - from: source node
#     - to: target node
#     - corr: signed correlation or similar strength in [-1, 1]
#
# Output:
#   A PDF, PNG, or SVG stress-layout network plot.
#
# Dependencies:
#   ggplot2, readr, ggraph, tidygraph, ggprism, graphlayouts, ggrepel
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (column names and labels).
#   Edit DATA PREPARATION to keep directed multi-edges or drop the
#   undirected collapse. Edit PLOT only when geometry or styling must change.
#
# Scientific assumptions:
#   Layout is computed from unweighted graph distance (stress). |corr| is
#   a visual edge width, not a layout weight (igraph treats weight as
#   distance, which would push strong links apart).
#   Reciprocal A-B / B-A rows are collapsed to the signed value with
#   larger |corr| so the layout is a simple undirected graph.
#   Colour encodes sign of corr. This script does not compute correlation
#   and does not run community detection.
#   Use graph-force / graph-basic / graph-cartesian when x/y are supplied.

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
        title = "Stress network",
        x = "",
        y = "",
        edge = "Sign",
        width = "|corr|"
    )
)

load_packages(c(
    "ggplot2", "readr", "ggraph", "tidygraph", "ggprism",
    "graphlayouts", "ggrepel"
))

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
self <- edges[[from_col]] == edges[[to_col]]
if (any(self)) {
    stop("Self-loops are not drawn by this template.", call. = FALSE)
}
edges$corr_abs <- abs(edges[[corr_col]])
# For each undirected pair, keep the signed row with larger |corr|.
pair_key <- paste(
    pmin(edges[[from_col]], edges[[to_col]]),
    pmax(edges[[from_col]], edges[[to_col]]),
    sep = "\t"
)
ord <- order(pair_key, -edges$corr_abs, seq_len(nrow(edges)))
edges <- edges[ord, , drop = FALSE]
pair_key <- pair_key[ord]
edges <- edges[!duplicated(pair_key), , drop = FALSE]
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
    directed = FALSE
)
set.seed(1)
layout <- ggraph::create_layout(graph, layout = "stress")

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
node_colour <- "#3369e7"
p <- ggraph::ggraph(layout) +
    ggraph::geom_edge_link(
        ggplot2::aes(colour = corr_label, width = corr_abs),
        alpha = 0.9
    ) +
    ggraph::geom_node_point(size = 3.8, colour = node_colour) +
    ggrepel::geom_text_repel(
        data = layout,
        ggplot2::aes(x = x, y = y, label = name),
        size = 4,
        colour = node_colour,
        inherit.aes = FALSE,
        seed = 1,
        max.overlaps = Inf,
        box.padding = 0.28,
        point.padding = 0.35,
        min.segment.length = 0,
        segment.size = 0.3,
        segment.colour = node_colour,
        segment.alpha = 0.55,
        force = 1.2,
        force_pull = 0.9
    ) +
    ggraph::scale_edge_colour_manual(
        name = config$labels$edge,
        values = c(
            Negative = "#ff4f81",
            Neutral = "#ffc168",
            Positive = "#2dde98"
        )
    ) +
    ggraph::scale_edge_width(
        name = config$labels$width,
        range = c(0.4, 2.4),
        limits = c(0, 1)
    ) +
    ggplot2::scale_x_continuous(expand = expansion(mult = 0.16)) +
    ggplot2::scale_y_continuous(expand = expansion(mult = 0.16)) +
    ggplot2::coord_fixed(clip = "off") +
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
        plot.margin = margin(12, 12, 12, 12)
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output, width = 8, height = 8)
