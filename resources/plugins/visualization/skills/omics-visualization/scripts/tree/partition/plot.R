#!/usr/bin/env Rscript

# Template-ID: tree-partition
#
# Purpose:
#   Draw a rectangular partition / icicle plot: each depth is a band,
#   and band width encodes leaf size.
#
# Inputs:
#   A hierarchical table (tsv/csv). Default example:
#     - parent: parent node name (empty for the root)
#     - name: node label
#     - path: unique node id
#     - value: leaf size (for example gene count)
#
# Output:
#   A PDF, PNG, or SVG icicle plot.
#
# Dependencies:
#   ggplot2, readr, ggraph, tidygraph, ggprism
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (column names and labels).
#   Edit DATA PREPARATION to change grouping, weights, or path parsing.
#   Edit PLOT only when the chart geometry or styling must change.
#
# Scientific assumptions:
#   path is the unique node id; parent names may repeat.
#   Band width encodes supplied leaf values. Internal nodes are layout
#   containers; their value is ignored and descendant leaves are summed
#   by the layout.
#   Values are treated as already summarized, not as raw replicates.
#   Fill is one catalog color per leaf term. Width still encodes the
#   supplied value. Depth is a layout coordinate, not a statistical
#   distance. Label size is fixed; it does not encode width.
#   This is the rectangular icicle. Use sunburst-basic for the circular
#   partition of the same table.

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
        title = "GO icicle",
        x = "",
        y = "",
        fill = "Term"
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
    stringsAsFactors = FALSE
)
has_child <- nodes$name %in% df$parent_id[!is.na(df$parent_id)]
nodes$is_leaf <- !has_child
nodes$weight <- ifelse(
    has_child,
    0,
    ifelse(is.na(nodes$value) | nodes$value <= 0, 1, nodes$value)
)
nodes$fill_id <- factor(
    nodes$label,
    levels = unique(nodes$label[nodes$is_leaf])
)
wrap_label <- function(x, width = 11) {
    vapply(
        as.character(x),
        function(s) paste(strwrap(s, width = width), collapse = "\n"),
        character(1),
        USE.NAMES = FALSE
    )
}
nodes$label <- ifelse(
    nodes$is_leaf,
    nodes$label,
    wrap_label(nodes$label, 11)
)
edges <- data.frame(
    from = df$parent_id[!is.na(df$parent_id)],
    to = df[[id_col]][!is.na(df$parent_id)],
    stringsAsFactors = FALSE
)
keep <- edges$from %in% nodes$name & edges$to %in% nodes$name
edges <- edges[keep, , drop = FALSE]
graph <- tidygraph::tbl_graph(nodes = nodes, edges = edges, directed = TRUE)
layout <- ggraph::create_layout(
    graph,
    layout = "partition",
    circular = FALSE,
    weight = weight
)

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
leaf_levels <- levels(droplevels(layout$fill_id[as.logical(layout$leaf)]))
fill_values <- palette_colors("Qualitative.Bold", n = length(leaf_levels))
label_values <- ifelse(.hex_luminance(fill_values) > 0.45, "#333333", "#ffffff")
layout$label_colour <- ifelse(
    as.logical(layout$leaf),
    unname(stats::setNames(label_values, leaf_levels)[as.character(layout$fill_id)]),
    "#333333"
)
layout$lab_angle <- ifelse(as.logical(layout$leaf), 90, 0)
y_top <- min(layout$y[as.integer(layout$depth) == 1L] - layout$height[as.integer(layout$depth) == 1L] / 2)
y_bottom <- max(layout$y[as.logical(layout$leaf)] + layout$height[as.logical(layout$leaf)] / 2)

p <- ggraph::ggraph(layout) +
    ggraph::geom_node_tile(
        ggplot2::aes(filter = !leaf & depth > 0),
        fill = "#e6e8ee",
        colour = "white",
        linewidth = 0.6
    ) +
    ggraph::geom_node_tile(
        ggplot2::aes(fill = fill_id, filter = leaf),
        colour = "white",
        linewidth = 0.6
    ) +
    ggraph::geom_node_text(
        ggplot2::aes(
            label = label,
            colour = label_colour,
            angle = lab_angle,
            filter = depth > 0
        ),
        size = 5.5,
        lineheight = 0.92
    ) +
    ggplot2::scale_fill_manual(
        name = config$labels$fill,
        values = stats::setNames(fill_values, leaf_levels),
        breaks = leaf_levels,
        na.value = "#e6e8ee"
    ) +
    ggplot2::scale_colour_identity(guide = "none") +
    ggplot2::scale_y_reverse() +
    ggplot2::coord_cartesian(ylim = c(y_bottom, y_top), clip = "off") +
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
        plot.margin = margin(12, 12, 12, 12)
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output, width = 12, height = 8)
