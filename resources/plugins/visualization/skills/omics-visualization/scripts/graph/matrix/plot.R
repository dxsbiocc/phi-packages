#!/usr/bin/env Rscript

# Template-ID: graph-matrix
#
# Purpose:
#   Draw a graph as an adjacency matrix of points. Rows and columns are
#   nodes; each edge is a point at their intersection.
#
# Inputs:
#   Two tables (tsv/csv):
#     1. nodes: name, group
#     2. links: from, to, corr
#
# Output:
#   A PDF, PNG, or SVG matrix plot.
#
# Dependencies:
#   ggplot2, readr, ggraph, tidygraph, ggprism
#
# Example:
#   Rscript plot.R nodes.tsv links.tsv output.pdf
#
# Agent adaptation:
#   For a new dataset, edit only CONFIG (column names and labels).
#   Pass the node table then the edge table on the CLI.
#   Group membership must be supplied; this script does not cluster nodes.
#   Edit DATA PREPARATION to change node order (default is group then name).
#   Edit PLOT only when the chart geometry or styling must change.
#
# Scientific assumptions:
#   group is a supplied grouping (for example pathway role), not a detected
#   community.
#   Node order is a display choice, not a clustering or spectral result.
#   corr is a supplied signed strength. Colour encodes sign; point size
#   encodes absolute magnitude. This script does not compute correlation.
#   The matrix is directed: (i, j) is from -> to. It is not mirrored.

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
        id = "name",
        group = "group"
    ),
    edge_columns = list(
        from = "from",
        to = "to",
        corr = "corr"
    ),
    labels = list(
        title = "Adjacency matrix",
        x = "",
        y = "",
        group = "Group",
        edge = "Sign",
        size = "|corr|"
    )
)

load_packages(c("ggplot2", "readr", "ggraph", "tidygraph", "ggprism"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
nodes <- read_table_auto(io$input[[1]])
require_columns(nodes, config$columns)
id_col <- config$columns$id
group_col <- config$columns$group
nodes[[id_col]] <- as.character(nodes[[id_col]])
nodes$node_group <- factor(
    as.character(nodes[[group_col]]),
    levels = unique(as.character(nodes[[group_col]]))
)
nodes$sort_key <- as.numeric(nodes$node_group) +
    rank(nodes[[id_col]], ties.method = "first") / (nrow(nodes) + 1)

links <- read_table_auto(io$input[[2]])
require_columns(links, config$edge_columns)
from_col <- config$edge_columns$from
to_col <- config$edge_columns$to
corr_col <- config$edge_columns$corr
links[[from_col]] <- as.character(links[[from_col]])
links[[to_col]] <- as.character(links[[to_col]])
links[[corr_col]] <- as.numeric(links[[corr_col]])
if (any(!is.finite(links[[corr_col]]))) {
    stop("corr must be numeric and finite.", call. = FALSE)
}
links$corr_abs <- abs(links[[corr_col]])
links$corr_label <- factor(
    ifelse(
        links[[corr_col]] > 0,
        "Positive",
        ifelse(links[[corr_col]] < 0, "Negative", "Neutral")
    ),
    levels = c("Negative", "Neutral", "Positive")
)

missing_nodes <- setdiff(
    unique(c(links[[from_col]], links[[to_col]])),
    nodes[[id_col]]
)
if (length(missing_nodes) > 0) {
    stop(
        paste(
            "Edge endpoints missing from the node table:",
            paste(missing_nodes, collapse = ", ")
        ),
        call. = FALSE
    )
}

graph <- tidygraph::tbl_graph(
    nodes = nodes,
    edges = data.frame(
        from = links[[from_col]],
        to = links[[to_col]],
        corr_label = links$corr_label,
        corr_abs = links$corr_abs,
        stringsAsFactors = FALSE
    ),
    node_key = id_col,
    directed = TRUE
)
layout <- ggraph::create_layout(
    graph,
    layout = "matrix",
    sort.by = sort_key
)

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
group_values <- c("#3be8b0", "#1aafd0", "#6a67ce", "#ffb900", "#fc636b")

label_offset <- min(layout$x, layout$y, na.rm = TRUE) - 0.85

p <- ggraph::ggraph(layout) +
    ggraph::geom_edge_point(
        ggplot2::aes(colour = corr_label, edge_size = corr_abs),
        mirror = FALSE
    ) +
    ggraph::geom_node_text(
        ggplot2::aes(
            x = label_offset,
            y = y,
            label = .data[[id_col]],
            colour = node_group
        ),
        hjust = 1,
        size = 4
    ) +
    ggraph::geom_node_text(
        ggplot2::aes(
            x = x,
            y = label_offset,
            label = .data[[id_col]],
            colour = node_group
        ),
        angle = 90,
        hjust = 1,
        vjust = 0.5,
        size = 4
    ) +
    ggplot2::scale_colour_manual(
        name = config$labels$group,
        values = expand_palette(group_values, nlevels(nodes$node_group))
    ) +
    ggraph::scale_edge_colour_manual(
        name = config$labels$edge,
        values = c(
            Negative = "#fc636b",
            Neutral = "#a8adb7",
            Positive = "#3be8b0"
        )
    ) +
    ggraph::scale_edge_size(
        name = config$labels$size,
        range = c(1.4, 7.5),
        limits = c(0, 1)
    ) +
    ggplot2::scale_x_continuous(expand = expansion(mult = 0.04, add = c(2.8, 0.6))) +
    ggplot2::scale_y_continuous(expand = expansion(mult = 0.04, add = c(2.8, 0.6))) +
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
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        plot.title = element_text(size = 16),
        plot.background = element_blank(),
        plot.margin = margin(12, 16, 16, 16)
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output, width = 9, height = 9)
