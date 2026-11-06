#!/usr/bin/env Rscript

# Template-ID: graph-module-panel
#
# Purpose:
#   Draw multi-panel co-occurrence networks side by side. Connected nodes
#   use a Fruchterman–Reingold layout; degree-zero isolates sit on an
#   outer ring. Node colour encodes supplied module; within-module edges
#   inherit that colour, cross-module edges are grey. Each panel caption
#   reports Nodes and Edges counts.
#
# Inputs:
#   Two tables (tsv/csv):
#     1. nodes: system, node, module  (optional kingdom)
#     2. links: system, from, to, edge_group  (optional sign, rho)
#
# Output:
#   A PDF, PNG, or SVG multi-panel module network figure.
#
# Dependencies:
#   ggplot2, readr, igraph, tidygraph, dplyr, tibble, patchwork, ggprism
#
# Example:
#   Rscript plot.R nodes.tsv links.tsv output.pdf
#
# Agent adaptation:
#   Edit only CONFIG (column names, palette, panel order, layout radii).
#   Pass the node table then the edge table on the CLI.
#   Edit PLOT only when the geometry must change.
#
# Scientific assumptions:
#   Module membership and edges are supplied. This script does not run
#   correlation, RMT thresholding, iDIRECT, community detection, or
#   topological metrics beyond vcount / ecount for captions.
#   FR is used only for connected nodes inside each panel; isolates are
#   placed on a polar ring. That hybrid layout is the figure claim —
#   it is not a general replacement for graph-stress.

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
        system = "system",
        node = "node",
        module = "module"
    ),
    edge_columns = list(
        system = "system",
        from = "from",
        to = "to",
        edge_group = "edge_group"
    ),
    # NULL → unique systems in data order.
    system_order = c("Mono.M", "Inter.M", "Inter.S", "Mono.S"),
    panel_letters = c("A", "B", "C", "D"),
    palettes = list(
        module = "Qualitative.Bold"
    ),
    # Named hex override; unmatched modules fall back to palette.
    module_colours = c(
        M1 = "#4C78A8",
        M2 = "#E45756",
        M3 = "#F1CE63",
        M4 = "#54A24B",
        M5 = "#2E5A88",
        M6 = "#B279A2",
        Isolate = "#D0D0D0",
        Inter = "#C8C8C8"
    ),
    isolate_module = "Isolate",
    inter_edge_group = "Inter",
    layout = list(
        fr_niter = 1500L,
        active_radius = 0.72,
        isolate_r_min = 0.88,
        isolate_r_max = 1.10,
        seed = 1L
    ),
    node_size = 1.35,
    edge_width = 0.28,
    edge_alpha = 0.72,
    isolate_colour = "#D0D0D0",
    isolate_size = 1.05
)

load_packages(c(
    "ggplot2", "readr", "igraph", "tidygraph", "dplyr", "tibble",
    "patchwork", "ggprism"
))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
nodes <- as.data.frame(read_table_auto(io$input[[1]]), stringsAsFactors = FALSE)
links <- as.data.frame(read_table_auto(io$input[[2]]), stringsAsFactors = FALSE)
require_columns(nodes, config$columns)
require_columns(links, config$edge_columns)

sys_col <- config$columns$system
node_col <- config$columns$node
mod_col <- config$columns$module
e_sys <- config$edge_columns$system
from_col <- config$edge_columns$from
to_col <- config$edge_columns$to
eg_col <- config$edge_columns$edge_group

nodes[[sys_col]] <- as.character(nodes[[sys_col]])
nodes[[node_col]] <- as.character(nodes[[node_col]])
nodes[[mod_col]] <- as.character(nodes[[mod_col]])
links[[e_sys]] <- as.character(links[[e_sys]])
links[[from_col]] <- as.character(links[[from_col]])
links[[to_col]] <- as.character(links[[to_col]])
links[[eg_col]] <- as.character(links[[eg_col]])

ok_n <- !is.na(nodes[[sys_col]]) & nzchar(nodes[[sys_col]]) &
    !is.na(nodes[[node_col]]) & nzchar(nodes[[node_col]])
nodes <- nodes[ok_n, , drop = FALSE]
ok_e <- !is.na(links[[e_sys]]) & nzchar(links[[e_sys]]) &
    !is.na(links[[from_col]]) & nzchar(links[[from_col]]) &
    !is.na(links[[to_col]]) & nzchar(links[[to_col]]) &
    links[[from_col]] != links[[to_col]]
links <- links[ok_e, , drop = FALSE]
if (!nrow(nodes)) {
    stop("No nodes to plot.", call. = FALSE)
}

systems <- config$system_order
if (is.null(systems)) {
    systems <- unique(nodes[[sys_col]])
} else {
    systems <- systems[systems %in% unique(nodes[[sys_col]])]
}
if (!length(systems)) {
    stop("No systems overlap between config$system_order and nodes.", call. = FALSE)
}
letters_use <- config$panel_letters
if (length(letters_use) < length(systems)) {
    letters_use <- c(letters_use, LETTERS)[seq_along(systems)]
}

# Module colour map across all panels.
mod_levels <- unique(c(
    nodes[[mod_col]][nodes[[mod_col]] != config$isolate_module],
    config$isolate_module,
    config$inter_edge_group
))
mod_levels <- mod_levels[!is.na(mod_levels) & nzchar(mod_levels)]
fixed <- config$module_colours
fixed <- fixed[names(fixed) %in% mod_levels]
need <- setdiff(mod_levels, names(fixed))
if (length(need)) {
    pal <- palette_colors(config$palettes$module, n = length(need))
    names(pal) <- need
    # Prefer grey for isolate / inter if not already fixed.
    if (config$isolate_module %in% need) {
        pal[[config$isolate_module]] <- config$isolate_colour
    }
    if (config$inter_edge_group %in% need) {
        pal[[config$inter_edge_group]] <- "#C8C8C8"
    }
    fixed <- c(fixed, pal)
}
mod_cols <- fixed[mod_levels]

layout_one <- function(nd, ed, seed) {
    # Node table first column must be the unique vertex name for igraph.
    vdf <- data.frame(
        name = nd[[node_col]],
        module = nd[[mod_col]],
        stringsAsFactors = FALSE
    )
    if (anyDuplicated(vdf$name)) {
        stop("Duplicate node names within a system.", call. = FALSE)
    }
    edf <- data.frame(
        from = ed[[from_col]],
        to = ed[[to_col]],
        edge_group = ed[[eg_col]],
        stringsAsFactors = FALSE
    )
    # Keep edges whose endpoints exist in this system's node table.
    keep <- edf$from %in% vdf$name & edf$to %in% vdf$name
    edf <- edf[keep, , drop = FALSE]

    g <- igraph::graph_from_data_frame(edf, directed = FALSE, vertices = vdf)
    deg <- igraph::degree(g)
    active <- names(deg)[deg > 0]
    iso <- names(deg)[deg == 0]

    xy <- matrix(0, nrow = igraph::vcount(g), ncol = 2)
    rownames(xy) <- igraph::V(g)$name

    set.seed(seed)
    if (length(active) >= 2L) {
        ga <- igraph::induced_subgraph(g, vids = active)
        la <- igraph::layout_with_fr(ga, niter = config$layout$fr_niter)
        # Centre and scale active cloud into the inner disc.
        la <- scale(la, center = TRUE, scale = FALSE)
        rmax <- max(sqrt(rowSums(la^2)), 1e-8)
        la <- la / rmax * config$layout$active_radius
        xy[active, ] <- la
    } else if (length(active) == 1L) {
        xy[active, ] <- c(0, 0)
    }

    if (length(iso)) {
        ang <- seq(0, 2 * pi, length.out = length(iso) + 1L)[-1]
        rr <- runif(
            length(iso),
            config$layout$isolate_r_min,
            config$layout$isolate_r_max
        )
        xy[iso, 1] <- rr * cos(ang)
        xy[iso, 2] <- rr * sin(ang)
    }

    list(graph = g, xy = xy, n_nodes = igraph::vcount(g), n_edges = igraph::ecount(g))
}

panel_plot <- function(sys, letter, seed) {
    nd <- nodes[nodes[[sys_col]] == sys, , drop = FALSE]
    ed <- links[links[[e_sys]] == sys, , drop = FALSE]
    # igraph wants Node as first vertex column conceptually; we build vdf above.
    lay <- layout_one(nd, ed, seed = seed)
    g <- lay$graph
    xy <- lay$xy

    tg <- tidygraph::as_tbl_graph(g)
    tg <- tg |>
        tidygraph::activate(nodes) |>
        dplyr::mutate(
            x = xy[name, 1],
            y = xy[name, 2],
            is_isolate = module == config$isolate_module | name %in% names(which(igraph::degree(g) == 0))
        )

    # Edge colour: named module → module colour; else inter grey.
    eg <- igraph::as_data_frame(g, what = "edges")
    eg$edge_colour <- ifelse(
        eg$edge_group %in% names(mod_cols) &
            !(eg$edge_group %in% c(config$isolate_module, config$inter_edge_group)),
        eg$edge_group,
        config$inter_edge_group
    )

    tg <- tg |>
        tidygraph::activate(edges) |>
        dplyr::mutate(edge_colour = eg$edge_colour)

    caption <- sprintf("Nodes:%d,  Edges:%d", lay$n_nodes, lay$n_edges)
    title <- sprintf("%s  %s", letter, sys)

    node_df <- tg |>
        tidygraph::activate(nodes) |>
        tibble::as_tibble()
    edge_df <- tg |>
        tidygraph::activate(edges) |>
        tibble::as_tibble()
    # Manual segments from layout coords (avoids re-layout).
    if (nrow(edge_df)) {
        edge_df$x <- node_df$x[edge_df$from]
        edge_df$y <- node_df$y[edge_df$from]
        edge_df$xend <- node_df$x[edge_df$to]
        edge_df$yend <- node_df$y[edge_df$to]
    }

    p <- ggplot2::ggplot() +
        ggplot2::geom_segment(
            data = edge_df,
            ggplot2::aes(
                x = x, y = y, xend = xend, yend = yend,
                colour = edge_colour
            ),
            linewidth = config$edge_width,
            alpha = config$edge_alpha,
            show.legend = FALSE
        ) +
        ggplot2::geom_point(
            data = node_df[node_df$is_isolate, , drop = FALSE],
            ggplot2::aes(x = x, y = y),
            colour = config$isolate_colour,
            size = config$isolate_size,
            stroke = 0
        ) +
        ggplot2::geom_point(
            data = node_df[!node_df$is_isolate, , drop = FALSE],
            ggplot2::aes(x = x, y = y, colour = module),
            size = config$node_size,
            stroke = 0,
            show.legend = FALSE
        ) +
        ggplot2::scale_colour_manual(values = mod_cols, guide = "none") +
        ggplot2::coord_fixed(xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2), clip = "off") +
        ggplot2::labs(title = title, caption = caption) +
        theme_prism(base_size = 10) +
        ggplot2::theme(
            axis.text = ggplot2::element_blank(),
            axis.ticks = ggplot2::element_blank(),
            axis.line = ggplot2::element_blank(),
            axis.title = ggplot2::element_blank(),
            panel.grid = ggplot2::element_blank(),
            plot.title = ggplot2::element_text(size = 11, face = "bold", hjust = 0),
            plot.caption = ggplot2::element_text(size = 9, hjust = 0.5, colour = "#333333"),
            plot.background = ggplot2::element_blank(),
            panel.background = ggplot2::element_blank(),
            plot.margin = ggplot2::margin(6, 4, 6, 4)
        )
    p
}

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
panels <- vector("list", length(systems))
for (i in seq_along(systems)) {
    panels[[i]] <- panel_plot(
        systems[[i]],
        letters_use[[i]],
        seed = config$layout$seed + i - 1L
    )
}
p <- patchwork::wrap_plots(panels, nrow = 1L)

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(
    p,
    io$output,
    width = max(10, 3.2 * length(systems)),
    height = 4.2
)
