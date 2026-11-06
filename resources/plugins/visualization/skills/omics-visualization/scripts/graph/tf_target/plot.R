#!/usr/bin/env Rscript

# Template-ID: graph-tf-target
#
# Purpose:
#   Draw a TF–target regulatory network: large labelled TF nodes, small
#   target nodes filled by supplied log2FC (diverging), edges coloured by
#   the source TF and transparent by supplied interaction weight
#   (preview.png; e.g. iRegulon / GENIE3 display).
#
# Inputs:
#   1. nodes.tsv: id, type (TF|Target), log2FC; optional x, y
#   2. links.tsv: source, target, weight
#
# Output:
#   A PDF, PNG, or SVG TF–target network plot.
#
# Dependencies:
#   ggplot2, readr, ggraph, tidygraph, ggprism, scales, ggrepel,
#   graphlayouts (when x/y are not supplied)
#
# Example:
#   Rscript plot.R nodes.tsv links.tsv output.pdf
#
# Agent adaptation:
#   Edit CONFIG for columns, sizes, fill limits, palettes, layout, and
#   labels. Extra TFs or targets are extra rows, not a new id. A plain
#   force hairball stays graph-force; a signed corr stress net stays
#   graph-stress. Prefer supplied x/y for paper-faithful layouts.
#
# Scientific assumptions:
#   Edges, weights, log2FC, and TF/Target labels are supplied. This
#   script does not run iRegulon, GENIE3, DE, or motif enrichment.
#   Layout without x/y is a display placement (stress), not a biological
#   claim.

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
        type = "type",
        log2fc = "log2FC",
        x = "x",
        y = "y"
    ),
    link_columns = list(
        source = "source",
        target = "target",
        weight = "weight"
    ),
    # Matched against type column (case-insensitive).
    tf_levels = c("TF", "Target"),
    size_values = c(TF = 11, Target = 3.2),
    fill_limits = c(-1, 1),
    palettes = list(
        fill = "Diverging.TealRose",
        edge = "Qualitative.Bold"
    ),
    edge_width = 0.35,
    edge_alpha_range = c(0.04, 0.55),
    # Used only when x/y are missing: "stress" | "fr".
    layout = "stress",
    layout_seed = 42,
    labels = list(
        title = NULL,
        fill = "log2FC",
        size = NULL,
        edge_alpha = "GENIE3 weight"
    ),
    size = list(width = 10, height = 9)
)

load_packages(c(
    "ggplot2", "readr", "ggraph", "tidygraph", "ggprism",
    "scales", "ggrepel", "graphlayouts"
))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
nodes <- as.data.frame(read_table_auto(io$input[[1]]), stringsAsFactors = FALSE)
links <- as.data.frame(read_table_auto(io$input[[2]]), stringsAsFactors = FALSE)

id_col <- config$columns$id
type_col <- config$columns$type
fc_col <- config$columns$log2fc
x_col <- config$columns$x
y_col <- config$columns$y
src_col <- config$link_columns$source
tgt_col <- config$link_columns$target
w_col <- config$link_columns$weight

require_columns(nodes, list(id = id_col, type = type_col, log2fc = fc_col))
require_columns(links, list(source = src_col, target = tgt_col, weight = w_col))

nodes[[id_col]] <- as.character(nodes[[id_col]])
nodes[[fc_col]] <- as.numeric(nodes[[fc_col]])
if (any(!is.finite(nodes[[fc_col]]))) {
    stop("log2FC must be numeric and finite.", call. = FALSE)
}
type_chr <- as.character(nodes[[type_col]])
type_norm <- ifelse(
    grepl("^tf$", type_chr, ignore.case = TRUE),
    "TF",
    "Target"
)
nodes$.type <- factor(type_norm, levels = c("TF", "Target"))
nodes$.size <- unname(config$size_values[as.character(nodes$.type)])

has_xy <- x_col %in% names(nodes) && y_col %in% names(nodes) &&
    all(is.finite(as.numeric(nodes[[x_col]]))) &&
    all(is.finite(as.numeric(nodes[[y_col]])))
if (has_xy) {
    nodes$.x <- as.numeric(nodes[[x_col]])
    nodes$.y <- as.numeric(nodes[[y_col]])
}

links[[src_col]] <- as.character(links[[src_col]])
links[[tgt_col]] <- as.character(links[[tgt_col]])
links[[w_col]] <- as.numeric(links[[w_col]])
links <- links[is.finite(links[[w_col]]) & links[[w_col]] > 0, , drop = FALSE]
if (!nrow(links)) stop("No finite positive-weight edges.", call. = FALSE)

node_ids <- nodes[[id_col]]
if (anyDuplicated(node_ids)) stop("Node id must be unique.", call. = FALSE)
miss <- setdiff(unique(c(links[[src_col]], links[[tgt_col]])), node_ids)
if (length(miss)) {
    stop(
        paste("Edge endpoints missing from nodes:", paste(head(miss, 8), collapse = ", ")),
        call. = FALSE
    )
}

# Edge colour key = source TF id.
tf_ids <- nodes[[id_col]][nodes$.type == "TF"]
links$.tf <- factor(links[[src_col]], levels = tf_ids)
if (any(is.na(links$.tf))) {
    stop("Every edge source must be a TF node.", call. = FALSE)
}
links$.weight <- links[[w_col]]

graph <- tidygraph::tbl_graph(
    nodes = nodes,
    edges = data.frame(
        from = links[[src_col]],
        to = links[[tgt_col]],
        .tf = links$.tf,
        .weight = links$.weight,
        stringsAsFactors = FALSE
    ),
    node_key = id_col,
    directed = TRUE
)

set.seed(config$layout_seed)
if (has_xy) {
    layout <- ggraph::create_layout(
        graph,
        layout = "manual",
        x = nodes$.x,
        y = nodes$.y
    )
} else {
    lay <- match.arg(config$layout, c("stress", "fr"))
    layout <- ggraph::create_layout(graph, layout = lay)
}

# Fill: palette runs cool (down) -> warm (up); mid = centre swatch.
fill_pal <- palette_colors(config$palettes$fill)
fill_low <- fill_pal[[1]]
fill_mid <- fill_pal[[ceiling(length(fill_pal) / 2)]]
fill_high <- fill_pal[[length(fill_pal)]]
fc_lim <- config$fill_limits

edge_cols <- expand_palette(
    palette_colors(config$palettes$edge),
    length(tf_ids)
)
names(edge_cols) <- tf_ids

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p <- ggraph::ggraph(layout) +
    ggraph::geom_edge_diagonal(
        ggplot2::aes(colour = .tf, alpha = .weight),
        width = config$edge_width,
        strength = 0.9,
        lineend = "round"
    ) +
    ggraph::scale_edge_colour_manual(
        values = edge_cols,
        guide = "none"
    ) +
    ggraph::scale_edge_alpha_continuous(
        name = config$labels$edge_alpha,
        range = config$edge_alpha_range,
        guide = ggplot2::guide_legend(
            title.position = "top",
            override.aes = list(colour = "grey35", width = 1.2),
            order = 3
        )
    ) +
    ggplot2::geom_point(
        data = layout,
        ggplot2::aes(
            x = x, y = y,
            fill = .data[[fc_col]],
            size = .type
        ),
        shape = 21,
        colour = "grey25",
        stroke = 0.25
    ) +
    ggplot2::scale_fill_gradient2(
        name = config$labels$fill,
        low = fill_low,
        mid = fill_mid,
        high = fill_high,
        midpoint = 0,
        limits = fc_lim,
        oob = scales::squish,
        guide = ggplot2::guide_colorbar(
            title.position = "top",
            barwidth = grid::unit(3.2, "cm"),
            barheight = grid::unit(0.35, "cm"),
            order = 1
        )
    ) +
    ggplot2::scale_size_manual(
        name = config$labels$size,
        values = config$size_values,
        breaks = c("TF", "Target"),
        labels = c("TFs", "Targets"),
        guide = ggplot2::guide_legend(
            title.position = "top",
            override.aes = list(
                fill = "white",
                colour = "grey25",
                stroke = 0.6
            ),
            order = 2
        )
    ) +
    ggrepel::geom_text_repel(
        data = layout[layout$.type == "TF", , drop = FALSE],
        ggplot2::aes(x = x, y = y, label = .data[[id_col]]),
        size = 3.2,
        fontface = "bold",
        colour = "grey10",
        box.padding = 0.35,
        point.padding = 0.25,
        min.segment.length = 0,
        segment.colour = "grey60",
        segment.size = 0.25,
        max.overlaps = Inf,
        seed = config$layout_seed
    ) +
    ggplot2::coord_equal(clip = "off") +
    ggplot2::labs(title = config$labels$title) +
    ggprism::theme_prism(base_size = 11) +
    ggplot2::theme(
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        axis.line = ggplot2::element_blank(),
        axis.title = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank(),
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.title = ggplot2::element_text(size = ggplot2::rel(0.85)),
        legend.text = ggplot2::element_text(size = ggplot2::rel(0.75)),
        plot.background = ggplot2::element_rect(fill = "white", colour = NA),
        panel.background = ggplot2::element_rect(fill = "white", colour = NA),
        plot.margin = ggplot2::margin(8, 8, 6, 8)
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(
    p,
    io$output,
    width = config$size$width,
    height = config$size$height
)
