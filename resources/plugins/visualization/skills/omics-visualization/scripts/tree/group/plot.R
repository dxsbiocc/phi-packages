#!/usr/bin/env Rscript

# Template-ID: tree-group
#
# Purpose:
#   Draw a forest of trees (multiple roots) from parent/path columns.
#
# Inputs:
#   A hierarchical table (tsv/csv). Default example:
#     - parent: parent node name (empty for roots)
#     - name: node label
#     - path: unique node id
#     - value: optional size
#
# Output:
#   A PDF, PNG, or SVG grouped tree plot.
#
# Dependencies:
#   ggplot2, readr, ggraph, tidygraph, ggprism
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (column names and labels).
#   Edit DATA PREPARATION to change category order or reshape the table.
#   Edit PLOT only when the chart geometry or styling must change.
#
# Scientific assumptions:
#   Empty parent rows are roots of separate trees.
#   path is the unique node id.
#   Layout is a display of the supplied parent-child edges.

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
        parent = "parent",
        label = "name",
        id = "path",
        value = "value"
    ),
    labels = list(
        title = "Grouped Tree",
        x = "",
        y = ""
    )
)

load_packages(c("ggplot2", "readr", "ggraph", "tidygraph", "ggprism"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
require_columns(df, config$columns)
id_col <- config$columns$id
parent_col <- config$columns$parent
label_col <- config$columns$label
value_col <- config$columns$value
df[[id_col]] <- as.character(df[[id_col]])
df[[parent_col]] <- as.character(df[[parent_col]])
df[[label_col]] <- as.character(df[[label_col]])
df[[value_col]] <- as.numeric(df[[value_col]])
root_missing <- is.na(df[[parent_col]]) | !nzchar(df[[parent_col]])
name_suffix <- paste0("/", df[[label_col]])
from_name <- substr(
    df[[id_col]],
    1L,
    nchar(df[[id_col]]) - nchar(name_suffix)
)
df$parent_id <- ifelse(
    root_missing,
    NA_character_,
    ifelse(
        endsWith(df[[id_col]], name_suffix),
        from_name,
        sub("/[^/]+$", "", df[[id_col]])
    )
)
df$parent_id[df$parent_id == df[[id_col]]] <- NA_character_

nodes <- data.frame(
    name = df[[id_col]],
    label = df[[label_col]],
    value = df[[value_col]],
    depth = lengths(strsplit(df[[id_col]], "/", fixed = TRUE)),
    stringsAsFactors = FALSE
)
if ("colour" %in% names(config$columns)) {
    colour_col <- config$columns$colour
    if (colour_col %in% names(df)) {
        nodes$colour <- as.character(df[[colour_col]])
        nodes$colour[!nzchar(nodes$colour) | is.na(nodes$colour)] <- NA_character_
    }
}
has_child <- nodes$name %in% df$parent_id[!is.na(df$parent_id)]
nodes$weight <- ifelse(has_child, 0, ifelse(is.na(nodes$value) | nodes$value <= 0, 1, nodes$value))
edges <- data.frame(
    from = df$parent_id[!is.na(df$parent_id)],
    to = df[[id_col]][!is.na(df$parent_id)],
    stringsAsFactors = FALSE
)
keep <- edges$from %in% nodes$name & edges$to %in% nodes$name
edges <- edges[keep, , drop = FALSE]
graph <- tidygraph::tbl_graph(nodes = nodes, edges = edges, directed = TRUE)


# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p <- ggraph::ggraph(graph, layout = "dendrogram", circular = FALSE) +
    ggraph::geom_edge_diagonal(colour = "#037ef3") +
    ggraph::geom_node_point(size = 1.4, colour = "#f85a40") +
    ggraph::geom_node_text(
        ggplot2::aes(label = label, filter = depth <= 2),
        hjust = 0,
        size = 3,
        nudge_y = 0.05
    ) +
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
save_ggplot(p, io$output, width = 10, height = 8)
