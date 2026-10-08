#!/usr/bin/env Rscript

# Template-ID: sunburst-basic
#
# Purpose:
#   Draw a sunburst chart of a hierarchy using path, parent, and value columns.
#
# Inputs:
#   A hierarchical table (tsv/csv). Default example:
#     - parent: parent node name (empty for roots)
#     - name: node label
#     - path: unique node id
#     - value: leaf or node size
#
# Output:
#   A PDF, PNG, or SVG sunburst chart.
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
#   path is the unique node id; parent names may repeat.
#   Internal-node size is the sum of descendant leaves in the layout.
#   Values are treated as already summarized, not as raw replicates.

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
        title = "Sunburst Basic",
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
    path_depth = lengths(strsplit(df[[id_col]], "/", fixed = TRUE)),
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
nodes$family <- vapply(
    strsplit(nodes$name, "/", fixed = TRUE),
    function(parts) parts[[1L]],
    character(1)
)
family_levels <- unique(nodes$family)
family_base <- expand_palette(
    c("#d7689d",
    "#f3d337",
    "#a0c443",
    "#66b5ae",
    "#5da4dc",
    "#40769e"),
    length(family_levels)
)
names(family_base) <- family_levels
lighten_hex <- function(hex, p) {
    rgb <- grDevices::col2rgb(hex) / 255
    rgb <- rgb + (1 - rgb) * p
    grDevices::rgb(rgb[1, ], rgb[2, ], rgb[3, ])
}
nodes$fill <- vapply(seq_len(nrow(nodes)), function(i) {
    lighten_hex(
        family_base[[nodes$family[[i]]]],
        min(0.4, (nodes$path_depth[[i]] - 1) * 0.12)
    )
}, character(1))
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
    circular = TRUE,
    weight = weight
)
mid_a <- (layout$start + layout$end) / 2
mid_r <- (layout$r0 + layout$r) / 2
layout$x <- mid_r * sin(mid_a)
layout$y <- mid_r * cos(mid_a)
lab_deg <- atan2(layout$y, layout$x) * 180 / pi
layout$lab_angle <- ifelse(lab_deg > 90 | lab_deg < -90, lab_deg + 180, lab_deg)
fill_vals <- stats::setNames(layout$fill, layout$name)

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p <- ggraph::ggraph(layout) +
    ggraph::geom_node_arc_bar(
        ggplot2::aes(fill = name),
        colour = "white",
        linewidth = 0.2
    ) +
    ggplot2::scale_fill_manual(values = fill_vals, guide = "none") +
    ggraph::geom_node_text(
        ggplot2::aes(label = label, angle = lab_angle),
        size = 3,
        colour = "#333333"
    ) +
    labs(
        title = NULL,
        caption = config$labels$title,
        x = NULL,
        y = NULL
    ) +
    theme_prism() +
    theme(
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none",
        plot.background = element_blank(),
        plot.caption = element_text(
            hjust = 0.5,
            face = "bold",
            size = 16,
            colour = "#3c3c41"
        )
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output, width = 8, height = 8)
