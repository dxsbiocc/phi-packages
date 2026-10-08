#!/usr/bin/env Rscript

# Template-ID: graph-force
#
# Purpose:
#   Draw a force-directed network using supplied node coordinates and edges.
#
# Inputs:
#   Two tables (tsv/csv):
#     1. nodes: id, name, x, y, symbolSize, category
#     2. links: source, target
#
# Output:
#   A PDF, PNG, or SVG force-layout network plot.
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
#   Node x/y are supplied layout coordinates, not recomputed here.
#   symbolSize is a supplied display size.
#   category is a supplied grouping, not a community detection result.

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
        id = "id",
        label = "name",
        x = "x",
        y = "y",
        size = "symbolSize",
        colour = "category"
    ),
    labels = list(
        title = "Force Graph",
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
label_col <- config$columns$label
x_col <- config$columns$x
y_col <- config$columns$y
size_col <- config$columns$size
colour_col <- config$columns$colour
nodes[[id_col]] <- as.character(nodes[[id_col]])
nodes[[x_col]] <- as.numeric(nodes[[x_col]])
nodes[[y_col]] <- as.numeric(nodes[[y_col]])
nodes[[size_col]] <- as.numeric(nodes[[size_col]])
nodes[[colour_col]] <- factor(nodes[[colour_col]], levels = unique(nodes[[colour_col]]))

links <- read_table_auto(io$input[[2]])
require_columns(links, list(source = "source", target = "target"))
links$from <- as.character(links$source)
links$to <- as.character(links$target)

graph <- tidygraph::tbl_graph(
    nodes = nodes,
    edges = links[, c("from", "to")],
    node_key = id_col,
    directed = FALSE
)

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p <- ggraph::ggraph(
    graph,
    layout = "manual",
    x = nodes[[x_col]],
    y = nodes[[y_col]]
) +
    ggraph::geom_edge_link(colour = "#86878c", alpha = 0.25, width = 0.3) +
    ggraph::geom_node_point(ggplot2::aes(
        size = .data[[size_col]],
        colour = .data[[colour_col]]
    )) +
    ggplot2::scale_size_continuous(range = c(1, 8)) +
    scale_colour_manual(
        values = expand_palette(
            c("#3fbe95",
            "#785db0",
            "#505372",
            "#fb628b",
            "#ff994d",
            "#b6d634",
            "#5070dd",
            "#0ca8df",
            "#ffd10a"),
            nlevels(nodes[[colour_col]])
        )
    ) +
    labs(title = config$labels$title, colour = "Group", size = "Size") +
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
