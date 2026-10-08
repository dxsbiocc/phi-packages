#!/usr/bin/env Rscript

# Template-ID: scatter-diagonal
#
# Purpose:
#   Draw a diagonal scatter plot of two log-transformed expression axes,
#   with labels for genes farthest from equal expression.
#
# Inputs:
#   A table with one row per gene. Default example:
#     - WT: wild-type expression
#     - KO: knockout expression
#     - group: direction label (None / Up / Down)
#     - symbol: gene label
#
# Output:
#   A PDF, PNG, or SVG diagonal scatter plot.
#
# Dependencies:
#   ggplot2, readr, ggprism, ggrepel
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (column names and labels).
#   Edit DATA PREPARATION to change the log2 transform or how many
#   labels are drawn per direction.
#   Edit PLOT only when the chart geometry or styling must change.
#
# Scientific assumptions:
#   Each row is one gene with two comparable expression values.
#   log2(WT) and log2(KO) are display transforms of raw expression.
#   Labels use log2(WT/KO) as a distance from the diagonal, not a
#   formal differential-expression test.

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
        x = "WT",
        y = "KO",
        colour = "group",
        label = "symbol"
    ),
    labels = list(
        title = "Yap1-KO",
        x = "WT",
        y = "KO"
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
colour_col <- config$columns$colour
label_col <- config$columns$label

df <- df[stats::complete.cases(df[, c(x_col, y_col, colour_col, label_col)]), ]
df[[x_col]] <- as.numeric(df[[x_col]])
df[[y_col]] <- as.numeric(df[[y_col]])
df[[colour_col]] <- factor(df[[colour_col]], levels = unique(df[[colour_col]]))
df$log2_x <- log2(pmax(df[[x_col]], 1e-12))
df$log2_y <- log2(pmax(df[[y_col]], 1e-12))

use_data <- df[tolower(as.character(df[[colour_col]])) != "none", , drop = FALSE]
use_data$log2FC <- log2(pmax(use_data[[x_col]], 1e-12) / pmax(use_data[[y_col]], 1e-12))
text_data <- do.call(rbind, lapply(
    split(use_data, use_data$log2FC > 0),
    function(sub) {
        sub[head(order(abs(sub$log2FC), decreasing = TRUE), 10), , drop = FALSE]
    }
))

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
colour_values <- c(
    Down = "#a16928",
    None = "#edeac2",
    Up = "#2887a1"
)

p <- ggplot(df, aes(
    x = log2_x,
    y = log2_y,
    colour = .data[[config$columns$colour]]
)) +
    geom_point(size = 3, alpha = 0.5, show.legend = FALSE) +
    geom_label_repel(
        data = text_data,
        aes(label = .data[[config$columns$label]]),
        max.overlaps = 2000,
        show.legend = FALSE
    ) +
    annotate(
        "label",
        x = c(-5, 10),
        y = c(10, -5),
        label = c("Up", "Down"),
        size = 6,
        fill = c("#a16928", "#2887a1"),
        colour = "#525e61",
        alpha = 0.8
    ) +
    scale_colour_manual(values = colour_values, drop = FALSE) +
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
