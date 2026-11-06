#!/usr/bin/env Rscript

# Template-ID: tree-daylight
#
# Purpose:
#   Draw a supplied phylogenetic tree in an unrooted daylight layout.
#   Tip points encode a categorical group (colour) and a numeric size.
#
# Inputs:
#   Two files:
#     1. Newick tree (.nwk / .tree / .newick)
#     2. tip table (tsv/csv): label, group, size
#
# Output:
#   A PDF, PNG, or SVG daylight tree with tip bubbles.
#
# Dependencies:
#   ggplot2, readr, ape, ggtree, tidytree
#
# Example:
#   Rscript plot.R tree.nwk tips.tsv output.pdf
#
# Agent adaptation:
#   Edit only CONFIG for tip column names, palette, and size range.
#   Pass the Newick path first, then the tip table, on the CLI.
#   Edit PLOT only when the geometry must change.
#
# Scientific assumptions:
#   Topology and tip labels come from the Newick file. This script does
#   not build or root the tree, and it does not run RedeR / TreeAndLeaf
#   force relaxation. branch.length = "none" is a display choice so
#   spacing follows topology, not branch lengths. Tip group and size
#   are supplied; colour and size scales are display only.

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

io <- parse_io_args(n_input = 2L, input_names = c("tree", "tips"))

# -----------------------------------------------------------------------------
# CONFIG  (edit this block for a new dataset)
# -----------------------------------------------------------------------------
config <- list(
    tip_columns = list(
        label = "label",
        group = "group",
        size = "size"
    ),
    labels = list(
        group = "Group",
        size = "Size"
    ),
    palettes = list(
        group = "Qualitative.Bold"
    ),
    # Display size range for tip bubbles (mapped from supplied size).
    size_range = c(1.8, 7.5),
    edge_colour = "grey25",
    edge_width = 0.28,
    tip_alpha = 1
)

load_packages(c("ggplot2", "readr", "ape", "ggtree", "tidytree"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
tree_path <- io$input[[1]]
if (!file.exists(tree_path)) {
    stop(paste("Tree file not found:", tree_path), call. = FALSE)
}
phy <- tryCatch(
    ape::read.tree(tree_path),
    error = function(e) {
        stop(paste("Could not read Newick tree:", conditionMessage(e)), call. = FALSE)
    }
)
if (inherits(phy, "multiPhylo")) {
    phy <- phy[[1]]
}
if (!inherits(phy, "phylo")) {
    stop("Input tree must be a Newick phylo object.", call. = FALSE)
}
n_tip <- ape::Ntip(phy)
if (n_tip < 3L) {
    stop("Need at least 3 tips for a daylight tree.", call. = FALSE)
}
phy$tip.label <- as.character(phy$tip.label)
if (anyDuplicated(phy$tip.label)) {
    stop("Tip labels in the Newick tree must be unique.", call. = FALSE)
}

tips <- as.data.frame(read_table_auto(io$input[[2]]), stringsAsFactors = FALSE)
require_columns(tips, config$tip_columns)
lab_col <- config$tip_columns$label
grp_col <- config$tip_columns$group
sz_col <- config$tip_columns$size
tips[[lab_col]] <- as.character(tips[[lab_col]])
tips[[grp_col]] <- as.character(tips[[grp_col]])
tips[[sz_col]] <- as.numeric(tips[[sz_col]])
if (anyDuplicated(tips[[lab_col]])) {
    stop("Tip labels in the tip table must be unique.", call. = FALSE)
}
if (any(!is.finite(tips[[sz_col]]) | tips[[sz_col]] <= 0)) {
    stop("size must be finite and > 0.", call. = FALSE)
}
missing_tips <- setdiff(phy$tip.label, tips[[lab_col]])
if (length(missing_tips)) {
    stop(
        paste(
            "Tip table is missing Newick tip labels:",
            paste(utils::head(missing_tips, 8L), collapse = ", ")
        ),
        call. = FALSE
    )
}
extra_tips <- setdiff(tips[[lab_col]], phy$tip.label)
if (length(extra_tips)) {
    tips <- tips[tips[[lab_col]] %in% phy$tip.label, , drop = FALSE]
}
# Keep tip-table order for first-seen group levels; remap to tip order.
tips <- tips[match(phy$tip.label, tips[[lab_col]]), , drop = FALSE]
group_levels <- unique(tips[[grp_col]])
tips[[grp_col]] <- factor(tips[[grp_col]], levels = group_levels)

anno <- data.frame(
    label = tips[[lab_col]],
    group = tips[[grp_col]],
    size = tips[[sz_col]],
    stringsAsFactors = FALSE
)
anno$group <- factor(anno$group, levels = group_levels)

phy_ann <- tidytree::full_join(phy, anno, by = "label")
# ggtree::fortify supplies daylight coordinates. Drawing via ggtree()
# itself can fail under ggplot2 4.x; build the figure with ggplot2.
lay <- ggtree::fortify(phy_ann, layout = "daylight", branch.length = "none")
lay$group <- factor(lay$group, levels = group_levels)
lay$size <- as.numeric(lay$size)

parent_xy <- lay[, c("node", "x", "y")]
names(parent_xy) <- c("parent", "x_parent", "y_parent")
edges <- merge(
    lay[, c("node", "parent", "x", "y")],
    parent_xy,
    by = "parent",
    all.x = TRUE
)
edges <- edges[edges$parent != edges$node & is.finite(edges$x_parent), , drop = FALSE]
tip_df <- lay[!is.na(lay$isTip) & as.logical(lay$isTip), , drop = FALSE]
if (!nrow(tip_df)) {
    stop("No tip coordinates returned by daylight layout.", call. = FALSE)
}

pal <- setNames(
    palette_colors(config$palettes$group, n = length(group_levels)),
    group_levels
)
sz_range <- as.numeric(config$size_range)
if (length(sz_range) != 2L || any(!is.finite(sz_range)) || sz_range[1] >= sz_range[2]) {
    stop("CONFIG$size_range must be two increasing finite numbers.", call. = FALSE)
}

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p <- ggplot() +
    geom_segment(
        data = edges,
        aes(x = x_parent, y = y_parent, xend = x, yend = y),
        colour = config$edge_colour,
        linewidth = config$edge_width,
        lineend = "round"
    ) +
    geom_point(
        data = tip_df,
        aes(x = x, y = y, size = size),
        colour = "white",
        show.legend = FALSE
    ) +
    geom_point(
        data = tip_df,
        aes(x = x, y = y, colour = group, size = size),
        alpha = config$tip_alpha
    ) +
    scale_colour_manual(
        values = pal,
        name = config$labels$group,
        breaks = group_levels
    ) +
    scale_size_continuous(
        range = sz_range,
        name = config$labels$size
    ) +
    guides(
        colour = guide_legend(order = 1, override.aes = list(size = 3.2, alpha = 1)),
        size = guide_legend(order = 2, override.aes = list(colour = "grey35", alpha = 1))
    ) +
    coord_equal(clip = "off") +
    theme_void(base_size = 10) +
    theme(
        legend.position = "right",
        legend.title = element_text(size = 9, colour = "black"),
        legend.text = element_text(size = 8, colour = "black"),
        plot.background = element_rect(fill = "white", colour = NA),
        plot.margin = margin(10, 10, 10, 10)
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output, width = 8.5, height = 7.2)
