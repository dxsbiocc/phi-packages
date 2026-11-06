#!/usr/bin/env Rscript

# Template-ID: tree-basic
#
# Purpose:
#   Draw a hierarchical tree from parent/path columns.
#
# Inputs:
#   A hierarchical table (tsv/csv). Default example:
#     - parent: parent node name (empty for the root)
#     - name: node label
#     - path: unique node id
#     - value: optional size
#
# Output:
#   A PDF, PNG, or SVG tree plot.
#
# Dependencies:
#   ggplot2, readr, ggraph, tidygraph, ggprism
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (column names and labels).
#   layout = "circular" (example_regulon.tsv) is the radial variant:
#   hidden root, branch colours, depth shapes, leaf size by value and
#   radial leaf labels, e.g. regulatory axis -> TF -> target regulons.
#   A target under several TFs is one leaf per path. Do not run
#   DoRothEA / decoupleR or DE in the script.
#   Edit DATA PREPARATION to change category order or reshape the table.
#   Edit PLOT only when the chart geometry or styling must change.
#
# Scientific assumptions:
#   path is the unique node id; parent names may repeat.
#   Layout is a display of the supplied parent-child edges.
#   value is a display size when present, not a statistical weight.

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
        title = "Tree Basic",
        subtitle = NULL,
        x = "",
        y = ""
    ),
    # "vertical" (default) or "circular": radial dendrogram with the
    # root hidden, nodes coloured by their depth-1 branch, shape by
    # depth, leaf size by value, radial leaf labels and a key.
    layout = "vertical",
    circular = list(
        hide_root = TRUE,
        # "Family.Name:index" per depth-1 branch; unlisted = Qualitative.Bold
        branch_colors = NULL,
        depth_shapes = c(16, 17, 16),        # per depth 1, 2, 3+
        depth_sizes = c(7, 5),               # depth 1, 2 (leaves use value)
        leaf_size_range = c(1.6, 5.5),
        depth_names = c("Regulatory axis", "Transcription factor",
            "Upregulated target gene"),
        branch_title = "Encoding key",
        shape_title = "Node meaning",
        size_note = "Target node size:\nranked DE signal",
        leaf_label_size = 2.9,
        inner_label_size = 3.6
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
if (identical(config$layout, "circular")) {
    cf <- config$circular
    par_of <- stats::setNames(df$parent_id, df[[id_col]])
    depth_of <- function(id) {
        d <- 0L
        while (!is.na(par_of[[id]])) {
            id <- par_of[[id]]
            d <- d + 1L
        }
        d
    }
    branch_of <- function(id) {
        prev <- id
        while (!is.na(par_of[[id]])) {
            prev <- id
            id <- par_of[[id]]
        }
        prev
    }
    nodes$depth <- vapply(nodes$name, depth_of, integer(1))
    nodes$branch <- nodes$label[match(vapply(nodes$name, branch_of,
        character(1)), nodes$name)]
    nodes$branch[nodes$depth == 0L] <- NA_character_
    br_lv <- unique(nodes$branch[nodes$depth == 1L])
    fallback <- palette_colors("Qualitative.Bold", n = length(br_lv))
    br_col <- stats::setNames(vapply(seq_along(br_lv), function(i) {
        ref <- cf$branch_colors[[br_lv[[i]]]]
        if (is.null(ref)) return(fallback[[i]])
        parts <- strsplit(ref, ":", fixed = TRUE)[[1]]
        palette_colors(parts[[1]])[[as.integer(parts[[2]])]]
    }, character(1)), br_lv)
    max_d <- max(nodes$depth)
    nodes$kind <- factor(pmin(nodes$depth, 3L), levels = 1:3,
        labels = cf$depth_names)
    lv <- nodes$value[nodes$depth == max_d]
    lr <- range(lv, na.rm = TRUE)
    sz <- cf$leaf_size_range
    nodes$size <- ifelse(nodes$depth == 1L, cf$depth_sizes[[1]],
        ifelse(nodes$depth == 2L, cf$depth_sizes[[2]],
            if (diff(lr) > 0) sz[[1]] + (nodes$value - lr[[1]]) / diff(lr) *
                diff(sz) else mean(sz)))
    graph <- tidygraph::tbl_graph(nodes = nodes, edges = edges, directed = TRUE)
    lay <- ggraph::create_layout(graph, layout = "dendrogram", circular = TRUE)
    lay$ang <- atan2(lay$y, lay$x) * 180 / pi
    lay$right <- abs(lay$ang) <= 90
    lay$lab_ang <- ifelse(lay$right, lay$ang, lay$ang + 180)
    show <- if (isTRUE(cf$hide_root)) lay$depth > 0 else rep(TRUE, nrow(lay))
    pts <- lay[show, , drop = FALSE]
    leaves <- pts[pts$leaf, , drop = FALSE]
    inner <- pts[!pts$leaf, , drop = FALSE]
    shapes <- stats::setNames(cf$depth_shapes, cf$depth_names)
    p <- ggraph::ggraph(lay) +
        ggraph::geom_edge_diagonal(colour = "grey86", edge_width = 0.7) +
        geom_point(data = pts, aes(x = x, y = y, colour = branch,
            shape = kind, size = size)) +
        geom_text(data = leaves, aes(x = x * 1.07, y = y * 1.07, label = label,
            angle = lab_ang, hjust = ifelse(right, 0, 1), colour = branch),
            size = cf$leaf_label_size, show.legend = FALSE) +
        geom_text(data = inner, aes(x = x, y = y, label = label,
            colour = branch), size = cf$inner_label_size, fontface = "bold",
            vjust = 2.1, show.legend = FALSE) +
        scale_colour_manual(values = br_col, breaks = br_lv,
            name = cf$branch_title, na.translate = FALSE) +
        scale_shape_manual(values = shapes, name = cf$shape_title,
            drop = FALSE) +
        scale_size_identity() +
        guides(colour = guide_legend(order = 1, override.aes = list(size = 5,
            shape = 16)), shape = guide_legend(order = 2,
            override.aes = list(size = c(5, 4, 2.6), colour = "grey35"))) +
        coord_fixed(xlim = c(-1.3, 1.3), ylim = c(-1.3, 1.3), clip = "off") +
        annotate("text", x = 1.5, y = -0.42, label = cf$size_note,
            hjust = 0, vjust = 1, size = 3.2) +
        labs(title = config$labels$title, subtitle = config$labels$subtitle) +
        theme_void(base_size = 12) +
        theme(
            plot.title = element_text(face = "bold", size = 16),
            plot.subtitle = element_text(colour = "grey35"),
            legend.title = element_text(face = "bold"),
            legend.text = element_text(face = "bold"),
            plot.background = element_blank(),
            plot.margin = margin(10, 10, 10, 10)
        )
    save_ggplot(p, io$output, width = 11, height = 9.5)
    quit(save = "no", status = 0)
}

p <- ggraph::ggraph(graph, layout = "dendrogram", circular = FALSE) +
    ggraph::geom_edge_diagonal(colour = "#1aafd0", edge_alpha = 0.7, edge_width = 0.4) +
    ggraph::geom_node_point(
        shape = 21,
        size = 2.2,
        fill = "white",
        colour = "#fc636b",
        stroke = 0.8
    ) +
    ggraph::geom_node_text(
        ggplot2::aes(label = label, filter = !leaf),
        size = 3,
        colour = "#333333",
        vjust = -0.9,
        hjust = 0.5
    ) +
    ggraph::geom_node_text(
        ggplot2::aes(label = label, filter = leaf),
        size = 2.4,
        colour = "#333333",
        angle = -90,
        hjust = 0,
        vjust = 0.5,
        nudge_y = -0.02
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.32, 0.12))) +
    labs(title = config$labels$title, x = NULL, y = NULL) +
    theme_prism() +
    theme(
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        panel.grid = element_blank(),
        plot.background = element_blank(),
        plot.margin = margin(8, 8, 8, 8)
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output, width = 14, height = 10)
