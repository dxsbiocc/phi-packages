#!/usr/bin/env Rscript

# Template-ID: graph-circular-concentric
#
# Purpose:
#   Draw a dual-ring concentric network: hub nodes on an inner circle and
#   category-grouped nodes on an outer circle, with weighted coloured edges.
#
# Inputs:
#   Two tables (tsv/csv):
#     1. nodes: id, label, ring, category, size
#        ring must be "inner" or "outer"
#     2. links: source, target, weight (positive integers)
#
# Output:
#   A PDF, PNG, or SVG concentric circular network plot.
#
# Dependencies:
#   ggplot2, readr, ggprism, ggforce
#
# Example:
#   Rscript plot.R nodes.tsv links.tsv output.pdf
#
# Agent adaptation:
#   For a new dataset, edit only CONFIG (column names, radii, labels, palette).
#   Pass the node table then the edge table on the CLI.
#   Edit PLOT only when the chart geometry or styling must change.
#
# Scientific assumptions:
#   Dual-ring placement is a layout choice, not an inferred hierarchy.
#   ring and category are supplied groupings, not community detection.
#   size and weight are supplied display encodings, not computed here.
#   Edges are drawn between supplied endpoints; this script does not
#   infer or filter biological interactions.

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
        label = "label",
        ring = "ring",
        category = "category",
        size = "size"
    ),
    edge_columns = list(
        source = "source",
        target = "target",
        weight = "weight"
    ),
    layout = list(
        inner_radius = 0.42,
        outer_radius = 1.00,
        category_gap = 0.35,
        label_offset = 0.14,
        category_label_offset = 0.32
    ),
    labels = list(
        title = "",
        node_class = "Node class",
        line_width = "Line width"
    ),
    palettes = list(
        class = "Qualitative.Bold"
    ),
    # Prefer cross-ring edges (outer ↔ inner). Same-ring edges are allowed
    # but coloured by the source node category.
    prefer_cross_ring = TRUE
)

load_packages(c("ggplot2", "readr", "ggprism", "ggforce"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
nodes <- as.data.frame(read_table_auto(io$input[[1]]), stringsAsFactors = FALSE)
require_columns(nodes, config$columns)
id_col <- config$columns$id
label_col <- config$columns$label
ring_col <- config$columns$ring
category_col <- config$columns$category
size_col <- config$columns$size

nodes[[id_col]] <- as.character(nodes[[id_col]])
nodes[[label_col]] <- as.character(nodes[[label_col]])
nodes[[ring_col]] <- tolower(as.character(nodes[[ring_col]]))
nodes[[category_col]] <- as.character(nodes[[category_col]])
nodes[[size_col]] <- as.numeric(nodes[[size_col]])

if (anyDuplicated(nodes[[id_col]])) {
    stop("Node id values must be unique.", call. = FALSE)
}
if (!all(nodes[[ring_col]] %in% c("inner", "outer"))) {
    stop('ring must be "inner" or "outer" for every node.', call. = FALSE)
}
if (any(!is.finite(nodes[[size_col]]) | nodes[[size_col]] <= 0)) {
    stop("size must be positive and finite.", call. = FALSE)
}

inner_idx <- which(nodes[[ring_col]] == "inner")
outer_idx <- which(nodes[[ring_col]] == "outer")
if (!length(inner_idx) || !length(outer_idx)) {
    stop("Need at least one inner node and one outer node.", call. = FALSE)
}

# Category order: hub / inner classes first (legend), then outer classes.
outer_cats <- unique(nodes[[category_col]][outer_idx])
inner_cats <- unique(nodes[[category_col]][inner_idx])
cat_levels <- unique(c(inner_cats, outer_cats))
nodes[[category_col]] <- factor(nodes[[category_col]], levels = cat_levels)

links <- as.data.frame(read_table_auto(io$input[[2]]), stringsAsFactors = FALSE)
require_columns(links, config$edge_columns)
src_col <- config$edge_columns$source
tgt_col <- config$edge_columns$target
weight_col <- config$edge_columns$weight
links[[src_col]] <- as.character(links[[src_col]])
links[[tgt_col]] <- as.character(links[[tgt_col]])
links[[weight_col]] <- as.numeric(links[[weight_col]])
if (any(!is.finite(links[[weight_col]]) | links[[weight_col]] <= 0)) {
    stop("weight must be positive and finite.", call. = FALSE)
}
if (any(abs(links[[weight_col]] - round(links[[weight_col]])) > 1e-8)) {
    stop(
        "weight must be positive integers for the discrete line-width legend.",
        call. = FALSE
    )
}
links[[weight_col]] <- as.integer(round(links[[weight_col]]))

missing_ends <- setdiff(
    unique(c(links[[src_col]], links[[tgt_col]])),
    nodes[[id_col]]
)
if (length(missing_ends) > 0) {
    stop(
        paste(
            "Edge endpoints missing from the node table:",
            paste(missing_ends, collapse = ", ")
        ),
        call. = FALSE
    )
}

ring_of <- setNames(nodes[[ring_col]], nodes[[id_col]])
cat_of <- setNames(as.character(nodes[[category_col]]), nodes[[id_col]])

src_ring <- unname(ring_of[links[[src_col]]])
tgt_ring <- unname(ring_of[links[[tgt_col]]])
cross_ring <- src_ring != tgt_ring
if (isTRUE(config$prefer_cross_ring) && !all(cross_ring)) {
    n_same <- sum(!cross_ring)
    warning(
        sprintf(
            "%d same-ring edge(s) present; edge colour falls back to the source category.",
            n_same
        ),
        call. = FALSE
    )
}

edge_category <- character(nrow(links))
for (i in seq_len(nrow(links))) {
    s <- links[[src_col]][[i]]
    t <- links[[tgt_col]][[i]]
    if (identical(ring_of[[s]], "outer")) {
        edge_category[[i]] <- cat_of[[s]]
    } else if (identical(ring_of[[t]], "outer")) {
        edge_category[[i]] <- cat_of[[t]]
    } else {
        edge_category[[i]] <- cat_of[[s]]
    }
}
links$edge_category <- factor(edge_category, levels = cat_levels)
weight_levels <- as.character(sort(unique(links[[weight_col]])))
links$weight_factor <- factor(
    as.character(links[[weight_col]]),
    levels = weight_levels
)
# ---- angular layout ---------------------------------------------------------
r_in <- config$layout$inner_radius
r_out <- config$layout$outer_radius
gap <- config$layout$category_gap

place_group <- function(n, start, span) {
    if (n <= 0) {
        return(numeric(0))
    }
    if (n == 1) {
        return(start + span / 2)
    }
    start + seq(0, span, length.out = n)
}

nodes$theta <- NA_real_
nodes$r <- ifelse(nodes[[ring_col]] == "inner", r_in, r_out)

# Outer: pack categories with equal within-group spacing and inter-group gaps.
n_outer <- length(outer_idx)
n_gap <- max(length(outer_cats), 1L)
usable <- 2 * pi - gap * n_gap
if (usable <= 0) {
    stop("category_gap is too large for the number of outer classes.", call. = FALSE)
}
# Weight span by group size so denser classes get more arc.
outer_sizes <- vapply(
    outer_cats,
    function(cat) sum(nodes[[category_col]][outer_idx] == cat),
    integer(1)
)
spans <- usable * outer_sizes / sum(outer_sizes)
cursor <- pi / 2
category_mids <- setNames(numeric(length(outer_cats)), outer_cats)
for (j in seq_along(outer_cats)) {
    cat <- outer_cats[[j]]
    idx <- which(
        nodes[[ring_col]] == "outer" &
            as.character(nodes[[category_col]]) == cat
    )
    span_j <- spans[[j]]
    thetas <- place_group(length(idx), cursor, span_j)
    nodes$theta[idx] <- thetas
    category_mids[[cat]] <- mean(range(thetas))
    cursor <- cursor + span_j + gap
}

# Inner: equal spacing, starting near the first outer category.
n_inner <- length(inner_idx)
inner_start <- pi / 2
nodes$theta[inner_idx] <- place_group(
    n_inner,
    inner_start,
    2 * pi * (n_inner - 1) / n_inner
)

nodes$x <- nodes$r * cos(nodes$theta)
nodes$y <- nodes$r * sin(nodes$theta)
nodes$is_inner <- nodes[[ring_col]] == "inner"

# Gene labels: outer sit slightly outside; inner sit just outside the hub.
label_r <- ifelse(
    nodes$is_inner,
    r_in + config$layout$label_offset * 0.85,
    r_out + config$layout$label_offset
)
nodes$lx <- label_r * cos(nodes$theta)
nodes$ly <- label_r * sin(nodes$theta)
# Keep text upright around the circle.
nodes$label_rot <- (nodes$theta * 180 / pi) %% 360
flip <- nodes$label_rot > 90 & nodes$label_rot < 270
nodes$label_rot[flip] <- nodes$label_rot[flip] + 180
nodes$label_hjust <- ifelse(flip, 1, 0)
# Category titles further out at group mid-angles.
cat_df <- data.frame(
    category = outer_cats,
    theta = as.numeric(category_mids),
    stringsAsFactors = FALSE
)
cat_r <- r_out + config$layout$category_label_offset
cat_df$x <- cat_r * cos(cat_df$theta)
cat_df$y <- cat_r * sin(cat_df$theta)
cat_df$rot <- (cat_df$theta * 180 / pi) %% 360
cat_flip <- cat_df$rot > 90 & cat_df$rot < 270
cat_df$rot[cat_flip] <- cat_df$rot[cat_flip] + 180
cat_df$hjust <- ifelse(cat_flip, 1, 0)

# Edge endpoints from node coordinates.
id_lookup <- match(links[[src_col]], nodes[[id_col]])
id_lookup_end <- match(links[[tgt_col]], nodes[[id_col]])
links$x <- nodes$x[id_lookup]
links$y <- nodes$y[id_lookup]
links$xend <- nodes$x[id_lookup_end]
links$yend <- nodes$y[id_lookup_end]
ring_guides <- data.frame(
    r = c(r_in, r_out),
    stringsAsFactors = FALSE
)

pal <- palette_colors(config$palettes$class, n = length(cat_levels))
names(pal) <- cat_levels

weight_levels <- levels(links$weight_factor)
weight_widths <- setNames(
    seq(0.35, 1.8, length.out = max(length(weight_levels), 1L)),
    weight_levels
)

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p <- ggplot() +
    geom_segment(
        data = links,
        aes(
            x = x,
            y = y,
            xend = xend,
            yend = yend,
            colour = edge_category,
            linewidth = weight_factor
        ),
        alpha = 0.72,
        lineend = "round"
    ) +
    ggforce::geom_circle(
        data = ring_guides,
        aes(x0 = 0, y0 = 0, r = r),
        inherit.aes = FALSE,
        colour = "#9aa0a8",
        linetype = "dashed",
        linewidth = 0.55
    ) +
    # Outer hollow nodes
    geom_point(
        data = nodes[!nodes$is_inner, , drop = FALSE],
        aes(
            x = x,
            y = y,
            colour = .data[[category_col]],
            size = .data[[size_col]]
        ),
        shape = 21,
        fill = "white",
        stroke = 1.35
    ) +
    # Inner filled hubs
    geom_point(
        data = nodes[nodes$is_inner, , drop = FALSE],
        aes(
            x = x,
            y = y,
            colour = .data[[category_col]],
            fill = .data[[category_col]],
            size = .data[[size_col]]
        ),
        shape = 21,
        stroke = 0.4
    ) +
    geom_text(
        data = nodes,
        aes(
            x = lx,
            y = ly,
            label = .data[[label_col]],
            colour = .data[[category_col]],
            angle = label_rot,
            hjust = label_hjust
        ),
        size = 2.7,
        vjust = 0.5,
        show.legend = FALSE
    ) +
    geom_text(
        data = cat_df,
        aes(x = x, y = y, label = category, angle = rot, hjust = hjust),
        colour = "#2b2f36",
        size = 3.4,
        fontface = "bold",
        vjust = 0.5,
        inherit.aes = FALSE
    ) +
    scale_colour_manual(
        name = config$labels$node_class,
        values = pal,
        breaks = cat_levels,
        drop = FALSE
    ) +
    scale_fill_manual(
        name = config$labels$node_class,
        values = pal,
        breaks = cat_levels,
        drop = FALSE,
        guide = "none"
    ) +
    scale_linewidth_manual(
        name = config$labels$line_width,
        values = weight_widths,
        breaks = weight_levels,
        drop = FALSE
    ) +
    scale_size_continuous(range = c(2.2, 9.5), guide = "none") +
    coord_fixed(clip = "off", xlim = c(-1.65, 1.65), ylim = c(-1.65, 1.65)) +
    labs(title = config$labels$title) +
    theme_prism(base_size = 12) +
    theme(
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "white", colour = NA),
        legend.position = "right",
        legend.title = element_text(size = 11, colour = "#1f2328"),
        legend.text = element_text(size = 10, colour = "#1f2328"),
        legend.background = element_rect(fill = "white", colour = NA),
        legend.key = element_rect(fill = "white", colour = NA),
        legend.box = "vertical",
        legend.box.margin = margin(0, 0, 0, 8),
        plot.title = element_text(size = 14),
        plot.background = element_rect(fill = "white", colour = NA),
        plot.margin = margin(16, 28, 16, 16)
    ) +
    guides(
        colour = guide_legend(
            order = 1,
            override.aes = list(
                shape = ifelse(cat_levels %in% inner_cats, 19, 21),
                fill = ifelse(
                    cat_levels %in% inner_cats,
                    unname(pal[cat_levels]),
                    "white"
                ),
                size = 3.8,
                stroke = 1.2,
                linewidth = 0,
                linetype = 0,
                alpha = 1
            )
        ),
        linewidth = guide_legend(order = 2)
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output, width = 11, height = 9)