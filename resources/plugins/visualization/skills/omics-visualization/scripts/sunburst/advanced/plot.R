#!/usr/bin/env Rscript

# Template-ID: sunburst-advanced
#
# Purpose:
#   Draw a sunburst chart using supplied per-node colours when present.
#
# Inputs:
#   A hierarchical table (tsv/csv). Default example:
#     - parent, name, path, value
#     - itemStyle.color: optional hex colour per node
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
#   path is the unique node id.
#   itemStyle.color is a supplied display colour, not computed here.
#   Internal-node size is the sum of descendant leaves in the layout.

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
        value = "value",
        colour = "itemStyle.color"
    ),
    labels = list(
        title = "",
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
parent_map <- stats::setNames(df$parent_id, df[[id_col]])
nodes$tree_level <- vapply(nodes$name, function(id) {
    n <- 1L
    seen <- character(0)
    cur <- id
    while (!is.na(parent_map[[cur]]) && nzchar(parent_map[[cur]])) {
        if (cur %in% seen) {
            break
        }
        seen <- c(seen, cur)
        cur <- parent_map[[cur]]
        n <- n + 1L
    }
    n
}, integer(1))
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
max_level <- max(layout$tree_level)
mid_n <- max(1L, max_level - 2L)
is_root <- layout$tree_level == 1L
is_leaf <- layout$tree_level == max_level
i_mid <- layout$tree_level - 1L
span <- 0.78 - 0.36
layout$r0 <- ifelse(
    is_root,
    0.16,
    ifelse(is_leaf, 0.78, 0.36 + (i_mid - 1) / mid_n * span)
)
layout$r <- ifelse(
    is_root,
    0.36,
    ifelse(is_leaf, 0.82, 0.36 + i_mid / mid_n * span)
)
mid_a <- (layout$start + layout$end) / 2
mid_r <- (layout$r0 + layout$r) / 2
r_text <- ifelse(is_leaf, layout$r + 0.045, mid_r)
layout$x <- r_text * sin(mid_a)
layout$y <- r_text * cos(mid_a)
rad_deg <- atan2(layout$y, layout$x) * 180 / pi
flip <- rad_deg > 90 | rad_deg < -90
layout$rad_angle <- ifelse(flip, rad_deg + 180, rad_deg)
tan_deg <- rad_deg - 90
tan_flip <- tan_deg > 90 | tan_deg < -90
layout$tan_angle <- ifelse(tan_flip, tan_deg + 180, tan_deg)
layout$out_hjust <- ifelse(flip, 1, 0)
fill_vals <- stats::setNames(layout$colour, layout$name)
fill_vals[is.na(fill_vals) | !nzchar(fill_vals)] <- "grey80"

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p <- ggraph::ggraph(layout) +
    ggraph::geom_node_arc_bar(
        ggplot2::aes(fill = name),
        colour = "white",
        linewidth = 0.25
    ) +
    ggplot2::scale_fill_manual(
        values = fill_vals, na.value = "grey80", guide = "none"
    ) +
    ggraph::geom_node_text(
        ggplot2::aes(
            label = label,
            angle = tan_angle,
            filter = tree_level == 1
        ),
        size = 3.2,
        colour = "#eeeeee"
    ) +
    ggraph::geom_node_text(
        ggplot2::aes(
            label = label,
            angle = rad_angle,
            filter = tree_level > 1 & tree_level < max_level
        ),
        size = 2.8,
        colour = "#333333"
    ) +
    ggraph::geom_node_text(
        ggplot2::aes(
            label = label,
            angle = rad_angle,
            hjust = out_hjust,
            filter = tree_level == max_level
        ),
        size = 2.6,
        colour = "#333333"
    ) +
    labs(title = NULL, x = NULL, y = NULL) +
    theme_prism() +
    theme(
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none",
        plot.background = element_blank(),
        plot.margin = margin(12, 12, 12, 12)
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output, width = 10, height = 10)
