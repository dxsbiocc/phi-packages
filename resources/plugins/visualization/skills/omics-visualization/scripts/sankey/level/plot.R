#!/usr/bin/env Rscript

# Template-ID: sankey-level
#
# Purpose:
#   Draw a multi-level Sankey diagram from a weighted edge list.
#
# Inputs:
#   A three-column table (tsv/csv). Default example:
#     - source: origin node
#     - target: destination node
#     - value: flow weight
#
# Output:
#   A PDF, PNG, or SVG multi-level Sankey diagram.
#
# Dependencies:
#   ggplot2, readr, ggsankey, ggprism
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
#   Each row is one directed flow with a supplied weight.
#   Node levels are the longest path from roots in the edge list.
#   Values are treated as already summarized, not as raw replicates.

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
        source = "source",
        target = "target",
        value = "value"
    ),
    labels = list(
        title = "Multi-level Sankey",
        x = "",
        y = ""
    )
)

load_packages(c("ggplot2", "readr", "ggsankey", "ggprism"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
require_columns(df, config$columns)
src_col <- config$columns$source
tgt_col <- config$columns$target
val_col <- config$columns$value
df[[src_col]] <- as.character(df[[src_col]])
df[[tgt_col]] <- as.character(df[[tgt_col]])
df[[val_col]] <- as.numeric(df[[val_col]])
df <- df[df[[val_col]] > 0, , drop = FALSE]

all_nodes <- unique(c(df[[src_col]], df[[tgt_col]]))
incoming <- split(df[[src_col]], df[[tgt_col]])
level <- setNames(rep(NA_integer_, length(all_nodes)), all_nodes)
roots <- setdiff(all_nodes, df[[tgt_col]])
level[roots] <- 0L
changed <- TRUE
guard <- 0L
while (changed && guard < 50L) {
    changed <- FALSE
    guard <- guard + 1L
    for (node in all_nodes) {
        parents <- incoming[[node]]
        if (is.null(parents)) next
        parent_levels <- level[parents]
        if (anyNA(parent_levels)) next
        new_level <- max(parent_levels) + 1L
        if (is.na(level[[node]]) || level[[node]] < new_level) {
            level[[node]] <- new_level
            changed <- TRUE
        }
    }
}
level[is.na(level)] <- 0L

plot_df <- data.frame(
    x = as.character(level[df[[src_col]]]),
    next_x = as.character(level[df[[tgt_col]]]),
    node = df[[src_col]],
    next_node = df[[tgt_col]],
    value = df[[val_col]],
    stringsAsFactors = FALSE
)
x_levels <- as.character(sort(unique(as.integer(c(plot_df$x, plot_df$next_x)))))
# Only sinks need terminal rows; copying every target double-counts
# intermediate nodes and makes bars taller than the attached flows.
sinks <- setdiff(all_nodes, unique(df[[src_col]]))
last <- plot_df[plot_df$next_node %in% sinks, , drop = FALSE]
last$x <- as.character(last$next_x)
last$node <- last$next_node
last$next_x <- NA_character_
last$next_node <- NA_character_
plot_df <- rbind(plot_df, last)
plot_df$x <- factor(plot_df$x, levels = x_levels)
plot_df$next_x <- factor(plot_df$next_x, levels = x_levels)
# Default space is dominated by the single-node root column; pad from root flow.
root_flow <- sum(df[[val_col]][df[[src_col]] %in% roots])
sankey_space <- root_flow * 0.008

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p <- ggplot(plot_df, aes(
    x = x,
    next_x = next_x,
    node = node,
    next_node = next_node,
    value = value,
    fill = node
)) +
    ggsankey::geom_sankey(
        flow.alpha = 0.5,
        node.color = "white",
        width = 0.15,
        space = sankey_space
    ) +
    scale_fill_manual(
        values = expand_palette(
            c("#a0c443",
            "#f3d337",
            "#66b5ae",
            "#d7689d"),
            length(unique(plot_df$node))
        )
    ) +
    labs(title = config$labels$title, fill = "") +
    theme_prism() +
    theme(
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        legend.position = "none",
        plot.background = element_blank()
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output, width = 10, height = 7)
