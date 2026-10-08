#!/usr/bin/env Rscript

# Template-ID: graph-basic
#
# Purpose:
#   Draw a directed node-link graph from supplied node coordinates and weighted edges.
#
# Inputs:
#   Two tables (tsv/csv):
#     1. nodes: name (id), x, y
#     2. links: source, target, value
#
# Output:
#   A PDF, PNG, or SVG network plot.
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
#   Edit PLOT only when the chart geometry or styling must change.
#
# Scientific assumptions:
#   Node x/y are supplied layout coordinates, not computed here.
#   Edge value is a display weight, not a statistical test.
#   Edges are treated as directed.

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
        x = "x",
        y = "y"
    ),
    labels = list(
        title = "Directed Graph",
        x = "",
        y = ""
    )
)

load_packages(c("ggplot2", "readr", "ggraph", "tidygraph", "ggprism"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
nodes <- read_table_auto(io$input[[1]])
require_columns(nodes, config$columns)
id_col <- config$columns$id
x_col <- config$columns$x
y_col <- config$columns$y
nodes[[x_col]] <- as.numeric(nodes[[x_col]])
nodes[[y_col]] <- as.numeric(nodes[[y_col]])
nodes[[id_col]] <- as.character(nodes[[id_col]])

links <- read_table_auto(io$input[[2]])
require_columns(links, list(source = "source", target = "target", value = "value"))
links$from <- as.character(links$source)
links$to <- as.character(links$target)
links$value <- as.numeric(links$value)

graph <- tidygraph::tbl_graph(
    nodes = nodes,
    edges = links[, c("from", "to", "value")],
    node_key = id_col,
    directed = TRUE
)

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p <- ggraph::ggraph(graph, layout = "manual", x = nodes[[x_col]], y = nodes[[y_col]]) +
    ggraph::geom_edge_link(
        ggplot2::aes(width = value),
        arrow = grid::arrow(length = grid::unit(4, "mm")),
        end_cap = ggraph::circle(4, "mm"),
        colour = "#86878c",
        alpha = 0.8
    ) +
    ggraph::geom_node_point(size = 8, colour = "#5070dd") +
    ggraph::geom_node_text(ggplot2::aes(label = .data[[id_col]]), colour = "white", size = 3) +
    ggraph::scale_edge_width(range = c(0.4, 1.8)) +
    labs(title = config$labels$title) +
    theme_prism() +
    theme(
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.grid = element_blank(),
        plot.background = element_blank()
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output, width = 8, height = 8)
