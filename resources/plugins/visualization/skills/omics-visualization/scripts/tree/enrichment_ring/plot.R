#!/usr/bin/env Rscript

# Template-ID: tree-enrichment-ring
#
# Purpose:
#   Draw a circular classification tree of enrichment terms, with a
#   quantitative track on the circumference (bars or points).
#   Outer geometry is config$outer, not a second template.
#
# Inputs:
#   One row per enriched term. Default example:
#     - from: supplied class (KEGG subclass or similar)
#     - to: term name
#     - Count: gene count from the enrichment table
#     - pvalue: supplied p or adjusted p
#
# Output:
#   A PDF, PNG, or SVG ring figure (circular tree inset in a polar track).
#
# Dependencies:
#   ggplot2, readr, ggraph, tidygraph, ggrepel, patchwork
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (column names, outer, hole, bar_width,
#   labels). Edit DATA PREPARATION to change class order. Edit PLOT only
#   when geometry must change. This figure is the source template: polar
#   geom_col/point plus ggraph circular dendrogram via inset_element.
#   Alignment is leaf order, not atan2: ggraph places leaves at equal
#   angles; the polar track must use the same order, the same clockwise
#   start at 12 o'clock, and expand = c(0, 0) so n bars occupy a full
#   360 degrees. Do not use discrete expand to open the 12 o'clock join:
#   that compresses every bar and walks them off the leaves. The 12
#   o'clock gap is the same (1 - bar_width) as between every other pair.
#   Do not cover it with a white rect: that paints a white fan on a
#   transparent canvas.
#
# Scientific assumptions:
#   from is a supplied classification, not a clustered or detected group.
#   Count and pvalue are supplied enrichment results; this script does
#   not run enrichment. -log10(pvalue) is a display transform.
#   Circular order follows the classification tree, not a statistical
#   ranking of terms. Angle does not encode Count (that is sunburst).
#   Use bar-enrichment-* for a ranked Cartesian list. Use sunburst when
#   ring angle should encode descendant size.

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
        from = "from",
        to = "to",
        count = "Count",
        pvalue = "pvalue"
    ),
    root = "KEGG",
    outer = "bar",
    inner_hole = 2.2,
    bar_width = 0.9,
    labels = list(
        title = "",
        fill = "Class",
        alpha = "-log10(p)"
    )
)

load_packages(c(
    "ggplot2", "readr", "ggraph", "tidygraph", "ggrepel", "patchwork"
))

wrap_label <- function(x, width = 15) {
    vapply(as.character(x), function(s) {
        paste(strwrap(s, width = width), collapse = "\n")
    }, character(1), USE.NAMES = FALSE)
}

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
require_columns(df, config$columns)
from_col <- config$columns$from
to_col <- config$columns$to
count_col <- config$columns$count
pvalue_col <- config$columns$pvalue
df[[from_col]] <- as.character(df[[from_col]])
df[[to_col]] <- as.character(df[[to_col]])
df[[count_col]] <- as.numeric(df[[count_col]])
df[[pvalue_col]] <- as.numeric(df[[pvalue_col]])
if (any(!nzchar(df[[from_col]]) | !nzchar(df[[to_col]]))) {
    stop("from and to must be non-empty.", call. = FALSE)
}
if (anyDuplicated(df[[to_col]])) {
    stop("Term names in the to column must be unique.", call. = FALSE)
}
if (any(!is.finite(df[[count_col]]) | df[[count_col]] <= 0)) {
    stop("Count must be positive and finite.", call. = FALSE)
}
if (any(!is.finite(df[[pvalue_col]]) | df[[pvalue_col]] < 0 | df[[pvalue_col]] > 1)) {
    stop("pvalue must be numeric in [0, 1].", call. = FALSE)
}
outer_geom <- match.arg(config$outer, c("bar", "point"))
bar_width <- config$bar_width
if (identical(outer_geom, "bar")) {
    if (!is.finite(bar_width) || bar_width <= 0 || bar_width > 1) {
        stop("config$bar_width must be in (0, 1].", call. = FALSE)
    }
}
root <- config$root
if (!nzchar(root)) {
    stop("config$root must be a non-empty node name.", call. = FALSE)
}
if (root %in% c(df[[from_col]], df[[to_col]])) {
    stop("config$root collides with a class or term name.", call. = FALSE)
}

# Same order as the source template: arrange(from, to), then vertices as
# root + unique(union(from, to)). ggraph walks that tree clockwise from
# 12 o'clock; the polar factor must keep this leaf sequence.
ord <- order(df[[from_col]], df[[to_col]])
df <- df[ord, , drop = FALSE]
subclass_levels <- unique(df[[from_col]])
node_names <- unique(c(df[[from_col]], df[[to_col]]))
vertices <- data.frame(
    name = c(root, node_names),
    level = c(
        0L,
        ifelse(node_names %in% subclass_levels, 1L, 2L)
    ),
    stringsAsFactors = FALSE
)
vertices$label_wrap <- ifelse(
    vertices$level == 1L,
    wrap_label(vertices$name, 15),
    vertices$name
)
edges <- rbind(
    data.frame(
        from = df[[from_col]],
        to = df[[to_col]],
        Count = df[[count_col]],
        pvalue = df[[pvalue_col]],
        stringsAsFactors = FALSE
    ),
    data.frame(
        from = root,
        to = subclass_levels,
        Count = 0,
        pvalue = 0,
        stringsAsFactors = FALSE
    )
)
edges$colour <- ifelse(edges$from == root, edges$to, edges$from)
edges$colour <- factor(edges$colour, levels = subclass_levels)
graph <- tidygraph::tbl_graph(
    vertices,
    edges,
    node_key = "name",
    directed = TRUE
)

leaf <- vertices[vertices$level == 2L, , drop = FALSE]
idx <- match(leaf$name, df[[to_col]])
if (anyNA(idx)) {
    stop("A tree leaf is missing from the enrichment table.", call. = FALSE)
}
n_leaf <- nrow(leaf)
leaf$from <- factor(df[[from_col]][idx], levels = subclass_levels)
leaf$Count <- df[[count_col]][idx]
leaf$pvalue <- df[[pvalue_col]][idx]
leaf$neglogp <- -log10(pmax(leaf$pvalue, 1e-300))
leaf$name <- factor(leaf$name, levels = leaf$name)
leaf$id <- seq_len(n_leaf)
# Bar centres sit at id; polar x runs 0.5 .. n+0.5, so the spoke is
# at (id - 0.5) / n, not id / n. Using id / n tilts every label.
theta <- 360 * (leaf$id - 0.5) / n_leaf
leaf$angle <- 90 - theta
leaf$hjust <- ifelse(leaf$angle < -90, 1, 0)
leaf$angle <- ifelse(leaf$angle < -90, leaf$angle + 180, leaf$angle)

class_colours <- palette_colors("Qualitative.Bold", length(subclass_levels))
names(class_colours) <- subclass_levels
hole <- max(leaf$Count) * config$inner_hole
label_pad <- max(leaf$Count) * 0.12

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
g <- ggraph::ggraph(graph, layout = "dendrogram", circular = TRUE) +
    ggraph::geom_edge_diagonal(
        ggplot2::aes(colour = colour),
        show.legend = FALSE
    ) +
    ggraph::geom_node_point(
        ggplot2::aes(filter = level == 1L, colour = name),
        size = 3,
        show.legend = FALSE
    ) +
    ggraph::geom_node_text(
        ggplot2::aes(
            filter = level == 1L,
            label = label_wrap,
            colour = name
        ),
        lineheight = 0.8,
        repel = TRUE,
        size = 3,
        show.legend = FALSE
    ) +
    ggplot2::coord_fixed(
        xlim = c(-1.2, 1.2),
        ylim = c(-1.2, 1.2),
        clip = "on"
    ) +
    ggraph::scale_edge_colour_manual(values = class_colours, guide = "none") +
    ggplot2::scale_colour_manual(values = class_colours, guide = "none") +
    ggplot2::labs(title = NULL) +
    theme_void() +
    theme(
        plot.background = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank(),
        plot.margin = margin(0, 0, 0, 0)
    )

if (identical(outer_geom, "bar")) {
    half <- bar_width / 2
    leaf$xmin <- leaf$id - half
    leaf$xmax <- leaf$id + half
    track <- ggplot2::ggplot(leaf, ggplot2::aes(x = id, y = Count)) +
        ggplot2::geom_rect(
            ggplot2::aes(
                xmin = xmin,
                xmax = xmax,
                ymin = 0,
                ymax = Count,
                fill = from,
                alpha = neglogp
            ),
            colour = NA,
            show.legend = FALSE
        ) +
        ggplot2::geom_text(
            ggplot2::aes(
                y = Count + label_pad,
                label = name,
                angle = angle,
                hjust = hjust
            ),
            size = 2,
            colour = "#333333"
        ) +
        ggplot2::scale_x_continuous(
            limits = c(0.5, n_leaf + 0.5),
            expand = c(0, 0)
        ) +
        ggplot2::scale_y_continuous(
            limits = c(-hole, NA),
            expand = c(0, 0)
        )
} else {
    track <- ggplot2::ggplot(leaf, ggplot2::aes(name, Count)) +
        ggplot2::geom_point(
            ggplot2::aes(
                y = 0,
                size = Count,
                fill = from,
                alpha = neglogp
            ),
            shape = 21,
            stroke = 0,
            show.legend = FALSE
        ) +
        ggplot2::geom_text(
            ggplot2::aes(
                y = 0.1,
                label = name,
                angle = angle,
                hjust = hjust
            ),
            size = 2,
            colour = "#333333"
        ) +
        ggplot2::scale_x_discrete(expand = c(0, 0)) +
        ggplot2::scale_y_continuous(
            limits = c(-config$inner_hole, 1),
            expand = c(0, 0)
        ) +
        ggplot2::scale_size_continuous(range = c(1, 4), guide = "none")
}

track <- track +
    ggplot2::coord_polar() +
    ggplot2::scale_fill_manual(values = class_colours, guide = "none") +
    ggplot2::scale_alpha_continuous(range = c(0.45, 1), guide = "none") +
    ggplot2::labs(title = NULL) +
    theme_void() +
    theme(
        plot.background = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank(),
        plot.margin = margin(0, 0, 0, 0),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank()
    )

p <- track + patchwork::inset_element(
    g,
    left = 0.2,
    bottom = 0.2,
    right = 0.8,
    top = 0.8,
    align_to = "full"
)
# patchwork paints an opaque white canvas unless the annotation
# theme fill is NA. ggplot2 4 themes warn here; the fill still applies.
p <- p + patchwork::plot_annotation(
    theme = ggplot2::theme(
        plot.background = ggplot2::element_rect(fill = NA, colour = NA)
    )
)

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
# ggplot2 4 S7 themes are not is_theme() for patchwork; fill = NA still works.
suppressWarnings(save_ggplot(p, io$output, width = 8, height = 8))
