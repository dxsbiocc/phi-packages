#!/usr/bin/env Rscript

# Template-ID: tree-dendrogram
#
# Purpose:
#   Cluster rows of a numeric matrix with dist + hclust and draw the tree.
#   Circular layout is config$circular, not a second template.
#
# Inputs:
#   A table with one leaf-id column, optional group, and numeric features.
#   Default example:
#     - sample: leaf identifier (one row per sample)
#     - group: supplied annotation (not a detected cluster)
#     - remaining columns: numeric features used for distance
#
# Output:
#   A PDF, PNG, or SVG clustering dendrogram.
#
# Dependencies:
#   ggplot2, readr, ggraph, tidygraph, ggprism
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (column names, distance, linkage,
#   circular, labels). Edit DATA PREPARATION to cluster columns instead
#   of rows, or to skip scaling. Edit PLOT only when geometry must change.
#
# Scientific assumptions:
#   This script does compute dist and hclust; that is the figure.
#   Default: z-score each feature (column), Euclidean distance, Ward.D2.
#   Branch height is the clustering height, not a p-value or bootstrap.
#   group is a supplied label painted onto leaves, not a cutree result.
#   Use tree-basic for a named parent/path hierarchy that is not clustered.

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
        id = "sample",
        group = "group"
    ),
    distance = "euclidean",
    method = "ward.D2",
    scale = TRUE,
    circular = FALSE,
    labels = list(
        title = "Sample dendrogram",
        group = "Group"
    )
)

load_packages(c("ggplot2", "readr", "ggraph", "tidygraph", "ggprism"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
require_columns(df, config$columns$id)
id_col <- config$columns$id
group_col <- config$columns$group
df[[id_col]] <- as.character(df[[id_col]])
if (anyDuplicated(df[[id_col]])) {
    stop("Leaf identifiers in the id column must be unique.", call. = FALSE)
}
skip <- c(id_col, group_col)
feature_cols <- setdiff(names(df), skip)
if (length(feature_cols) < 2L) {
    stop("Need at least two numeric feature columns for distance.", call. = FALSE)
}
if (nrow(df) < 2L) {
    stop("Need at least two rows (leaves) to cluster.", call. = FALSE)
}
for (col in feature_cols) {
    df[[col]] <- as.numeric(df[[col]])
}
mat <- as.matrix(as.data.frame(df[feature_cols], optional = TRUE))
rownames(mat) <- df[[id_col]]
if (any(!is.finite(mat))) {
    stop("Feature columns must be numeric and finite.", call. = FALSE)
}
if (isTRUE(config$scale)) {
    mat <- scale(mat, center = TRUE, scale = TRUE)
    if (any(!is.finite(mat))) {
        stop(
            "Scaling produced non-finite values (a feature has zero variance).",
            call. = FALSE
        )
    }
}
hc <- stats::hclust(
    stats::dist(mat, method = config$distance),
    method = config$method
)
circular <- isTRUE(config$circular)
layout <- ggraph::create_layout(
    hc,
    layout = "dendrogram",
    circular = circular,
    height = height
)
has_group <- !is.null(group_col) && group_col %in% names(df)
if (has_group) {
    group_levels <- unique(as.character(df[[group_col]]))
    group_map <- stats::setNames(as.character(df[[group_col]]), df[[id_col]])
    layout$group <- factor(
        unname(group_map[layout$label]),
        levels = group_levels
    )
} else {
    layout$group <- factor("leaf")
}

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
group_n <- nlevels(droplevels(layout$group[layout$leaf]))
group_values <- palette_colors("Qualitative.Safe", max(group_n, 1L))
height_max <- max(layout$height, na.rm = TRUE)
if (!is.finite(height_max) || height_max <= 0) {
    height_max <- 1
}

leaf_df <- as.data.frame(layout)
leaf_df <- leaf_df[leaf_df$leaf, , drop = FALSE]

p <- ggraph::ggraph(layout) +
    ggraph::geom_edge_elbow(
        colour = "#6699cc",
        edge_width = 0.55,
        edge_alpha = 0.9
    ) +
    ggplot2::geom_point(
        data = leaf_df,
        ggplot2::aes(x = x, y = y, colour = group),
        size = 3.2,
        inherit.aes = FALSE
    )

if (circular) {
    p <- p +
        ggplot2::geom_text(
            data = leaf_df,
            ggplot2::aes(
                x = x,
                y = y,
                label = label,
                colour = group,
                angle = ggraph::node_angle(x, y)
            ),
            hjust = "outward",
            size = 4,
            inherit.aes = FALSE,
            show.legend = FALSE
        ) +
        ggplot2::coord_fixed(clip = "off") +
        ggplot2::scale_x_continuous(expand = expansion(mult = 0.28)) +
        ggplot2::scale_y_continuous(expand = expansion(mult = 0.28))
} else {
    p <- p +
        ggplot2::geom_text(
            data = leaf_df,
            ggplot2::aes(x = x, y = y, label = label, colour = group),
            angle = 90,
            hjust = 1,
            nudge_y = -height_max * 0.03,
            size = 4,
            inherit.aes = FALSE,
            show.legend = FALSE
        ) +
        ggplot2::coord_cartesian(clip = "off") +
        ggplot2::scale_y_continuous(
            expand = expansion(mult = c(0.28, 0.08))
        )
}

if (has_group) {
    p <- p + ggplot2::scale_colour_manual(
        name = config$labels$group,
        values = group_values,
        na.translate = FALSE
    )
} else {
    p <- p + ggplot2::scale_colour_manual(
        values = c(leaf = "#333333"),
        guide = "none"
    )
}

p <- p +
    ggplot2::labs(
        title = config$labels$title,
        x = NULL,
        y = NULL
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
        plot.margin = if (circular) {
            margin(16, 16, 16, 16)
        } else {
            margin(12, 16, 72, 16)
        }
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
if (circular) {
    save_ggplot(p, io$output, width = 8, height = 8)
} else {
    save_ggplot(p, io$output, width = 10, height = 6)
}
