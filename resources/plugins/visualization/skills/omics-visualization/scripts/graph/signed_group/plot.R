#!/usr/bin/env Rscript

# Template-ID: graph-signed-group
#
# Purpose:
#   MetaNet-style signed multi-omics / cross-domain network. Connected
#   components are laid out locally then packed on a circle (bubble of
#   islands). Node fill = class, node shape = layer (e.g. circle vs
#   square for Bacteria vs Fungi). Edges: colour = corr sign, width =
#   |corr|, linetype = inter vs intra layer. Covers three-omics and
#   bacteria–fungi with the same geometry.
#
# Inputs:
#   Two tables (tsv/csv):
#     1. nodes: id, class, size; optional layer, label
#     2. links: from, to, corr in [-1, 1]
#
# Output:
#   A PDF, PNG, or SVG annotated signed network.
#
# Dependencies:
#   ggplot2, readr, igraph, ggraph, tidygraph, ggprism, ggrepel
#
# Example:
#   Rscript plot.R nodes.tsv links.tsv output.pdf
#   Rscript plot.R nodes_bf.tsv links_bf.tsv bacteria_fungi.png
#
# Agent adaptation:
#   Edit only CONFIG. Pass nodes then links on the CLI.
#
# Scientific assumptions:
#   Edges and annotations are supplied. No Spearman / MetaNet /
#   community detection. Component packing is display geometry, not a
#   new layout-algorithm template id. Bacteria–fungi vs three-omics is
#   CONFIG (tables), not a new id.

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
        class = "class",
        size = "size",
        layer = "layer",
        label = "id"
    ),
    edge_columns = list(
        from = "from",
        to = "to",
        corr = "corr"
    ),
    labels = list(
        title = "Multi-omics network",
        class = NULL,          # NULL → use layer name as legend title when one class/layer
        size = "Abundance",
        edge = "Correlation",
        width = "abs(r)",
        relation = "Omics relationship"
    ),
    palettes = list(
        class = "Qualitative.Safe"
    ),
    # Match MetaNet-style: positive blue-ish / negative red-ish as in many
    # published panels; override in CONFIG if the paper flips them.
    edge_colours = c(
        Positive = "#4C72B0",
        Negative = "#C44E52",
        Neutral = "#B0B0B0"
    ),
    # Shape per layer (21 circle, 22 square, 24 triangle…). Unlisted → 21.
    layer_shapes = c(
        Bacteria = 21,
        Microbiome = 21,
        Fungi = 22,
        fugi = 22,
        Metabolome = 22,
        Transcriptome = 24
    ),
    # Named hex for common MetaNet-style classes; others fall back to palette.
    class_colours = c(
        "D_0__Bacteria" = "#7CB82F",
        Bacteria = "#7CB82F",
        fugi = "#F0C419",
        Fungi = "#F0C419",
        "k__Fungi" = "#9B6B9E"
    ),
    labels_num = 8L,
    label_colour = "#C44E52",
    node_size_range = c(2.2, 8.5),
    edge_width_range = c(0.35, 1.8),
    edge_curvature = 0.15,
    # Pack connected components on a circle (MetaNet multi-island look).
    component_radius = 1.0,
    component_gap = 0.35,
    fr_niter = 500L,
    seed = 1L,
    legend_number = TRUE
)

load_packages(c(
    "ggplot2", "readr", "igraph", "ggraph", "tidygraph", "ggprism", "ggrepel"
))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
nodes <- as.data.frame(read_table_auto(io$input[[1]]), stringsAsFactors = FALSE)
links <- as.data.frame(read_table_auto(io$input[[2]]), stringsAsFactors = FALSE)
require_columns(nodes, config$columns[c("id", "class", "size")])
require_columns(links, config$edge_columns)

id_col <- config$columns$id
class_col <- config$columns$class
size_col <- config$columns$size
layer_col <- config$columns$layer
label_col <- config$columns$label
from_col <- config$edge_columns$from
to_col <- config$edge_columns$to
corr_col <- config$edge_columns$corr

nodes[[id_col]] <- as.character(nodes[[id_col]])
nodes[[class_col]] <- as.character(nodes[[class_col]])
nodes[[size_col]] <- as.numeric(nodes[[size_col]])
if (anyDuplicated(nodes[[id_col]])) {
    stop("Duplicate node ids.", call. = FALSE)
}
if (any(!is.finite(nodes[[size_col]]) | nodes[[size_col]] < 0)) {
    stop("size must be finite and non-negative.", call. = FALSE)
}
has_layer <- !is.null(layer_col) && layer_col %in% names(nodes)
if (has_layer) {
    nodes[[layer_col]] <- as.character(nodes[[layer_col]])
} else {
    layer_col <- ".layer"
    nodes[[layer_col]] <- "Group"
    has_layer <- TRUE
}
if (is.null(label_col) || !label_col %in% names(nodes)) {
    label_col <- id_col
}
nodes[[label_col]] <- as.character(nodes[[label_col]])

links[[from_col]] <- as.character(links[[from_col]])
links[[to_col]] <- as.character(links[[to_col]])
links[[corr_col]] <- as.numeric(links[[corr_col]])
ok <- is.finite(links[[corr_col]]) &
    links[[from_col]] != links[[to_col]] &
    links[[from_col]] %in% nodes[[id_col]] &
    links[[to_col]] %in% nodes[[id_col]]
links <- links[ok, , drop = FALSE]
if (!nrow(links)) {
    stop("No finite edges with known endpoints.", call. = FALSE)
}
links$corr_abs <- abs(links[[corr_col]])
pair_key <- paste(
    pmin(links[[from_col]], links[[to_col]]),
    pmax(links[[from_col]], links[[to_col]]),
    sep = "\t"
)
ord <- order(pair_key, -links$corr_abs, seq_len(nrow(links)))
links <- links[ord, , drop = FALSE]
pair_key <- pair_key[ord]
links <- links[!duplicated(pair_key), , drop = FALSE]

links$corr_label <- factor(
    ifelse(
        links[[corr_col]] > 0, "Positive",
        ifelse(links[[corr_col]] < 0, "Negative", "Neutral")
    ),
    levels = c("Positive", "Negative", "Neutral")
)

# inter = endpoints in different layers; intra = same layer.
layer_of <- setNames(nodes[[layer_col]], nodes[[id_col]])
links$relation <- ifelse(
    layer_of[links[[from_col]]] == layer_of[links[[to_col]]],
    "intra",
    "inter"
)
links$relation <- factor(links$relation, levels = c("inter", "intra"))

keep_ids <- unique(c(links[[from_col]], links[[to_col]]))
nodes <- nodes[nodes[[id_col]] %in% keep_ids, , drop = FALSE]
nodes[[class_col]] <- factor(nodes[[class_col]], levels = unique(nodes[[class_col]]))
nodes[[layer_col]] <- factor(nodes[[layer_col]], levels = unique(nodes[[layer_col]]))

# Shape from layer.
shape_map <- config$layer_shapes
nodes$`.shape` <- unname(vapply(as.character(nodes[[layer_col]]), function(ly) {
    if (!is.null(shape_map) && ly %in% names(shape_map)) {
        as.integer(shape_map[[ly]])
    } else {
        21L
    }
}, integer(1)))

# ---- component-packed layout (circular islands) ----
g <- igraph::graph_from_data_frame(
    data.frame(
        from = links[[from_col]],
        to = links[[to_col]],
        stringsAsFactors = FALSE
    ),
    directed = FALSE,
    vertices = data.frame(name = nodes[[id_col]], stringsAsFactors = FALSE)
)
comps <- igraph::components(g)
memb <- comps$membership
comp_ids <- sort(unique(memb))
n_comp <- length(comp_ids)

set.seed(as.integer(config$seed))
xy <- matrix(0, nrow = igraph::vcount(g), ncol = 2)
rownames(xy) <- igraph::V(g)$name

# Component centres on a circle; size scales with sqrt(n).
comp_sizes <- vapply(comp_ids, function(cid) sum(memb == cid), integer(1))
ang <- seq(0, 2 * pi, length.out = n_comp + 1L)[-1]
# Slight radius jitter so rings don't look perfectly regular.
r_pack <- as.numeric(config$component_radius) *
    (0.55 + 0.45 * sqrt(comp_sizes / max(comp_sizes)))
centres <- cbind(r_pack * cos(ang), r_pack * sin(ang))

for (i in seq_along(comp_ids)) {
    cid <- comp_ids[[i]]
    vids <- names(memb)[memb == cid]
    if (length(vids) == 1L) {
        xy[vids, ] <- centres[i, ]
        next
    }
    sub <- igraph::induced_subgraph(g, vids = vids)
    if (length(vids) == 2L) {
        la <- rbind(c(-0.15, 0), c(0.15, 0))
    } else {
        la <- igraph::layout_with_fr(sub, niter = as.integer(config$fr_niter))
        la <- scale(la, center = TRUE, scale = FALSE)
        rmax <- max(sqrt(rowSums(la^2)), 1e-8)
        # Local radius grows slowly with component size.
        local_r <- 0.12 + 0.08 * sqrt(length(vids)) +
            as.numeric(config$component_gap) * 0.15
        la <- la / rmax * local_r
    }
    xy[vids, ] <- sweep(la, 2, centres[i, ], `+`)
}

# Normalise whole canvas to roughly unit circle.
xy <- scale(xy, center = TRUE, scale = FALSE)
rmax <- max(sqrt(rowSums(xy^2)), 1e-8)
xy <- xy / rmax

nodes$x <- xy[nodes[[id_col]], 1]
nodes$y <- xy[nodes[[id_col]], 2]

# Edge endpoints for geom_segment / curve.
links$x <- nodes$x[match(links[[from_col]], nodes[[id_col]])]
links$y <- nodes$y[match(links[[from_col]], nodes[[id_col]])]
links$xend <- nodes$x[match(links[[to_col]], nodes[[id_col]])]
links$yend <- nodes$y[match(links[[to_col]], nodes[[id_col]])]

# Top-N labels by size.
n_lab <- as.integer(config$labels_num)
nodes$`.label` <- ""
if (is.finite(n_lab) && n_lab > 0L) {
    ord_sz <- order(nodes[[size_col]], decreasing = TRUE)
    take <- ord_sz[seq_len(min(n_lab, length(ord_sz)))]
    nodes$`.label`[take] <- nodes[[label_col]][take]
}

# Colours.
class_levels <- levels(nodes[[class_col]])
if (!is.null(config$class_colours) && length(config$class_colours)) {
    class_cols <- config$class_colours
    miss <- setdiff(class_levels, names(class_cols))
    if (length(miss)) {
        extra <- palette_colors(config$palettes$class, n = length(miss))
        names(extra) <- miss
        class_cols <- c(class_cols, extra)
    }
    class_cols <- class_cols[class_levels]
} else {
    class_cols <- palette_colors(config$palettes$class, n = length(class_levels))
    names(class_cols) <- class_levels
}

# Legend labels with optional counts.
legend_lab <- function(levels_chr, counts) {
    if (isTRUE(config$legend_number)) {
        sprintf("%s: %d", levels_chr, counts[levels_chr])
    } else {
        levels_chr
    }
}
class_counts <- table(nodes[[class_col]])
class_legend <- legend_lab(class_levels, class_counts)
names(class_legend) <- class_levels

corr_counts <- table(links$corr_label)
corr_levels <- levels(links$corr_label)
corr_levels <- corr_levels[corr_levels %in% names(corr_counts)]
corr_legend <- legend_lab(as.character(corr_levels), corr_counts)
names(corr_legend) <- corr_levels

rel_counts <- table(links$relation)
rel_levels <- levels(links$relation)
rel_legend <- legend_lab(as.character(rel_levels), rel_counts)
names(rel_legend) <- rel_levels

# Build a combined aesthetic key for fill+shape legend grouped by layer:
# MetaNet shows bacteria / fungi as separate legend blocks. We approximate
# with one fill legend (class) and one shape legend (layer).
layer_levels <- levels(nodes[[layer_col]])
shape_vals <- unique(nodes[, c(layer_col, ".shape"), drop = FALSE])
shape_map_plot <- setNames(shape_vals$`.shape`, as.character(shape_vals[[layer_col]]))

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
curv <- as.numeric(config$edge_curvature)
p <- ggplot2::ggplot() +
    ggplot2::geom_curve(
        data = links,
        ggplot2::aes(
            x = x, y = y, xend = xend, yend = yend,
            colour = corr_label,
            linewidth = corr_abs,
            linetype = relation
        ),
        curvature = curv,
        alpha = 0.7,
        lineend = "round"
    ) +
    ggplot2::geom_point(
        data = nodes,
        ggplot2::aes(
            x = x, y = y,
            size = .data[[size_col]],
            fill = .data[[class_col]],
            shape = .data[[layer_col]]
        ),
        colour = "grey25",
        stroke = 0.3
    ) +
    ggrepel::geom_text_repel(
        data = nodes[nzchar(nodes$`.label`), , drop = FALSE],
        ggplot2::aes(x = x, y = y, label = .data$`.label`),
        size = 2.6,
        colour = config$label_colour,
        inherit.aes = FALSE,
        seed = as.integer(config$seed),
        max.overlaps = Inf,
        box.padding = 0.2,
        point.padding = 0.15,
        min.segment.length = 0,
        segment.size = 0.2,
        segment.colour = "grey60"
    ) +
    ggplot2::scale_colour_manual(
        name = config$labels$edge,
        values = config$edge_colours[corr_levels],
        breaks = corr_levels,
        labels = unname(corr_legend[corr_levels]),
        guide = ggplot2::guide_legend(order = 2, override.aes = list(linewidth = 1.2, linetype = "solid"))
    ) +
    ggplot2::scale_linetype_manual(
        name = config$labels$relation,
        values = c(inter = "solid", intra = "dashed"),
        breaks = rel_levels,
        labels = unname(rel_legend[rel_levels]),
        guide = ggplot2::guide_legend(order = 4, override.aes = list(colour = "grey30", linewidth = 0.8))
    ) +
    ggplot2::scale_linewidth(
        name = config$labels$width,
        range = as.numeric(config$edge_width_range),
        guide = ggplot2::guide_legend(order = 3)
    ) +
    ggplot2::scale_fill_manual(
        name = if (!is.null(config$labels$class)) config$labels$class else "Class",
        values = class_cols,
        breaks = class_levels,
        labels = unname(class_legend[class_levels]),
        guide = ggplot2::guide_legend(
            order = 1,
            override.aes = list(
                shape = 21,
                size = 3.5,
                colour = "grey25"
            )
        )
    ) +
    ggplot2::scale_shape_manual(
        name = "Layer",
        values = shape_map_plot,
        guide = ggplot2::guide_legend(order = 5, override.aes = list(size = 3.5, fill = "grey70"))
    ) +
    ggplot2::scale_size_continuous(
        name = config$labels$size,
        range = as.numeric(config$node_size_range),
        guide = ggplot2::guide_legend(order = 6, override.aes = list(shape = 21, fill = "grey70"))
    ) +
    ggplot2::coord_fixed(clip = "off") +
    ggplot2::labs(title = config$labels$title) +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(
        plot.title = ggplot2::element_text(
            size = 14, face = "bold", hjust = 0.5, colour = "#222222",
            margin = ggplot2::margin(b = 8)
        ),
        legend.title = ggplot2::element_text(size = 9, face = "bold"),
        legend.text = ggplot2::element_text(size = 8),
        legend.key.height = grid::unit(0.35, "cm"),
        legend.spacing.y = grid::unit(0.15, "cm"),
        plot.background = ggplot2::element_rect(fill = "white", colour = NA),
        plot.margin = ggplot2::margin(12, 12, 12, 12)
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output, width = 12, height = 7.5)
