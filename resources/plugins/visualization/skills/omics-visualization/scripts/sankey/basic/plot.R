#!/usr/bin/env Rscript

# Template-ID: sankey-basic
#
# Purpose:
#   Draw a Sankey diagram of weighted flows from source to target.
#
# Inputs:
#   A three-column table (tsv/csv). Default example:
#     - source: origin node
#     - target: destination node
#     - value: flow weight
#
# Output:
#   A PDF, PNG, or SVG Sankey diagram.
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
#   Node positions are a layout of the edge list, not a clustering result.
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
        title = "Sankey Diagram",
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
plot_df <- rbind(
    data.frame(
        x = 1,
        next_x = 2,
        node = df[[src_col]],
        next_node = df[[tgt_col]],
        value = df[[val_col]],
        stringsAsFactors = FALSE
    ),
    data.frame(
        x = 2,
        next_x = NA_real_,
        node = df[[tgt_col]],
        next_node = NA_character_,
        value = df[[val_col]],
        stringsAsFactors = FALSE
    )
)

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
    ggsankey::geom_sankey(flow.alpha = 0.6, node.color = "white") +
    scale_fill_manual(
        values = expand_palette(
            c("#e15f5f",
            "#944995",
            "#5468ce",
            "#2f9d9e",
            "#099d84"),
            length(unique(plot_df$node))
        )
    ) +
    ggsankey::geom_sankey_text(aes(label = node), size = 3, color = "black") +
    scale_x_continuous(breaks = c(1, 2), labels = c("source", "target")) +
    labs(title = config$labels$title, fill = "") +
    theme_prism() +
    theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        plot.background = element_blank()
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output, width = 8, height = 6)
