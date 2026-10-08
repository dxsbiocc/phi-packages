#!/usr/bin/env Rscript

# Template-ID: graph-hive
#
# Purpose:
#   Draw a hive plot: nodes sit on categorical radial axes, edges curve
#   between axes.
#
# Inputs:
#   Two tables (tsv/csv):
#     1. nodes: name, axis
#     2. links: from, to, corr
#
# Output:
#   A PDF, PNG, or SVG hive plot.
#
# Dependencies:
#   ggplot2, readr, ggraph, tidygraph, ggprism, ggrepel
#
# Example:
#   Rscript plot.R nodes.tsv links.tsv output.pdf
#
# Agent adaptation:
#   For a new dataset, edit only CONFIG (column names and labels).
#   Pass the node table then the edge table on the CLI.
#   Axis membership must be supplied; this script does not cluster nodes.
#   Edit PLOT only when the chart geometry or styling must change.
#
# Scientific assumptions:
#   axis is a supplied grouping (for example pathway role), not a detected
#   community.
#   corr is a supplied signed strength. Colour encodes sign; edge width
#   encodes absolute magnitude. This script does not compute correlation.
#   Hive placement is a layout choice, not a spatial embedding.

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
        axis = "axis"
    ),
    edge_columns = list(
        from = "from",
        to = "to",
        corr = "corr"
    ),
    labels = list(
        title = "Hive plot",
        x = "",
        y = "",
        axis = "Axis",
        edge = "Sign",
        width = "|corr|"
    )
)

load_packages(c("ggplot2", "readr", "ggraph", "tidygraph", "ggprism", "ggrepel"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
nodes <- read_table_auto(io$input[[1]])
require_columns(nodes, config$columns)
id_col <- config$columns$id
axis_col <- config$columns$axis
nodes[[id_col]] <- as.character(nodes[[id_col]])
nodes$hive_axis <- factor(
    as.character(nodes[[axis_col]]),
    levels = unique(as.character(nodes[[axis_col]]))
)

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
layout <- ggraph::create_layout(graph, layout = "hive", axis = hive_axis)

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
axis_values <- c("#3be8b0", "#1aafd0", "#6a67ce", "#ffb900", "#fc636b")

p <- ggraph::ggraph(layout) +
    ggraph::geom_axis_hive(
        ggplot2::aes(colour = hive_axis),
        linewidth = 1.2,
        label = FALSE
    ) +
    ggraph::geom_edge_hive(
        ggplot2::aes(colour = corr_label, width = corr_abs),
        alpha = 0.88
    ) +
    ggraph::geom_node_point(
        ggplot2::aes(colour = hive_axis),
        size = 2.8
    ) +
    ggrepel::geom_text_repel(
        data = layout,
        ggplot2::aes(
            x = x,
            y = y,
            label = .data[[id_col]],
            colour = hive_axis
        ),
        size = 3,
        inherit.aes = FALSE,
        seed = 42,
        max.overlaps = Inf,
        box.padding = 0.28,
        point.padding = 0.35,
        min.segment.length = 0,
        segment.size = 0.3,
        segment.alpha = 0.6,
        force = 1.4,
        force_pull = 0.8
    ) +
    ggplot2::scale_colour_manual(
        name = config$labels$axis,
        values = expand_palette(axis_values, nlevels(nodes$hive_axis))
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
        range = c(0.35, 2.2),
        limits = c(0, 1)
    ) +
    ggplot2::scale_x_continuous(expand = expansion(mult = 0.22)) +
    ggplot2::scale_y_continuous(expand = expansion(mult = 0.22)) +
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
        plot.margin = margin(18, 28, 18, 18)
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output, width = 9, height = 9)
