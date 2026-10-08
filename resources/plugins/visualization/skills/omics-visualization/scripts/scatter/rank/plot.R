#!/usr/bin/env Rscript

# Template-ID: scatter-rank
#
# Purpose:
#   Draw a rank scatter plot of log2 fold-change against gene rank, with
#   labels at both extremes of the ranking.
#
# Inputs:
#   A table with one row per gene. Default example:
#     - rank: rank order
#     - log2FC: log2 fold-change
#     - symbol: gene label
#
# Output:
#   A PDF, PNG, or SVG rank scatter plot.
#
# Dependencies:
#   ggplot2, readr, ggprism, ggrepel
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (column names and labels).
#   Edit DATA PREPARATION to change how many genes are labelled at
#   each end of the rank axis.
#   Edit PLOT only when the chart geometry or styling must change.
#
# Scientific assumptions:
#   Each row is one gene with a rank and a log2 fold-change.
#   Rank is treated as already computed; this script does not re-rank.
#   Labels mark the lowest and highest ranks, not a statistical cutoff.

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
        x = "rank",
        y = "log2FC",
        label = "symbol"
    ),
    labels = list(
        title = "Yap1-KO",
        x = "Rank",
        y = "log2FC"
    )
)

load_packages(c("ggplot2", "readr", "ggprism", "ggrepel"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
require_columns(df, config$columns)

x_col <- config$columns$x
y_col <- config$columns$y
label_col <- config$columns$label

df <- df[stats::complete.cases(df[, c(x_col, y_col, label_col)]), ]
df[[x_col]] <- as.numeric(df[[x_col]])
df[[y_col]] <- as.numeric(df[[y_col]])

topk <- 10
text_data <- rbind(
    df[head(order(df[[x_col]], decreasing = TRUE), topk), , drop = FALSE],
    df[head(order(df[[x_col]], decreasing = FALSE), topk), , drop = FALSE]
)

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p <- ggplot(df, aes(
    x = .data[[config$columns$x]],
    y = .data[[config$columns$y]],
    colour = .data[[config$columns$y]]
)) +
    geom_point(
        aes(size = abs(.data[[config$columns$y]])),
        stroke = 0,
        show.legend = FALSE
    ) +
    geom_point(
        data = text_data,
        size = 5,
        fill = "transparent",
        colour = "grey60",
        shape = 21,
        show.legend = FALSE
    ) +
    geom_label_repel(
        data = text_data,
        aes(label = .data[[config$columns$label]]),
        max.overlaps = 2000,
        show.legend = FALSE
    ) +
    scale_colour_gradient2(low = "#7fbc41", mid = "#f7f7f7", high = "#de77ae") +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y
    ) +
    theme_prism() +
    theme(
        plot.background = element_blank(),
        legend.title = element_text(size = 12)
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output)
