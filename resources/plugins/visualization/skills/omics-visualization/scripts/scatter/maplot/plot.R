#!/usr/bin/env Rscript

# Template-ID: scatter-maplot
#
# Purpose:
#   Draw an MA plot of average expression versus log2 fold-change, with
#   labels for the strongest up- and down-regulated genes.
#
# Inputs:
#   A table with one row per gene. Default example:
#     - baseMean: mean expression
#     - log2FC: log2 fold-change
#     - group: direction label (None / Up / Down)
#     - symbol: gene label
#
# Output:
#   A PDF, PNG, or SVG MA plot.
#
# Dependencies:
#   ggplot2, readr, ggprism, ggrepel
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (column names and labels).
#   Edit DATA PREPARATION to change the log10 transform or how many
#   labels are drawn per direction.
#   Edit PLOT only when the chart geometry or styling must change.
#
# Scientific assumptions:
#   Each row is one gene with a mean expression and a log2 fold-change.
#   log10(baseMean) is a display transform of average expression.
#   Labels are taken from non-"None" genes with the largest |log2FC|.

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
        x = "baseMean",
        y = "log2FC",
        group = "group",
        label = "symbol"
    ),
    labels = list(
        title = "Yap1-KO",
        x = "log10(baseMean)",
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
group_col <- config$columns$group
label_col <- config$columns$label

df <- df[stats::complete.cases(df[, c(x_col, y_col, group_col, label_col)]), ]
df[[x_col]] <- as.numeric(df[[x_col]])
df[[y_col]] <- as.numeric(df[[y_col]])
df[[group_col]] <- factor(df[[group_col]], levels = unique(df[[group_col]]))
df$log10_mean <- log10(pmax(df[[x_col]], 1e-12))

use_data <- df[tolower(as.character(df[[group_col]])) != "none", , drop = FALSE]
text_data <- do.call(rbind, lapply(
    split(use_data, use_data[[y_col]] < 0),
    function(sub) {
        sub[head(order(abs(sub[[y_col]]), decreasing = TRUE), 10), , drop = FALSE]
    }
))

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
colour_values <- c(Down = "#a3ad62", Up = "#df91a3")
fill_values <- c(Down = "#a3ad62", None = "grey70", Up = "#df91a3")
shape_values <- c(Down = 25, None = 21, Up = 24)

p <- ggplot(df, aes(
    x = log10_mean,
    y = .data[[config$columns$y]],
    label = .data[[config$columns$label]]
)) +
    geom_point(
        aes(
            fill = .data[[config$columns$group]],
            shape = .data[[config$columns$group]]
        ),
        alpha = 0.7,
        size = 3,
        colour = "grey50",
        stroke = 0.1,
        show.legend = FALSE
    ) +
    geom_point(
        data = text_data,
        aes(colour = .data[[config$columns$group]]),
        size = 5,
        fill = "transparent",
        shape = 21,
        show.legend = FALSE
    ) +
    geom_label_repel(
        data = text_data,
        aes(colour = .data[[config$columns$group]]),
        size = 4,
        max.overlaps = 2000,
        show.legend = FALSE
    ) +
    annotate(
        "label",
        x = c(5, 5),
        y = c(6, -6),
        label = c("Up", "Down"),
        size = 6,
        fill = c("#df91a3", "#a3ad62"),
        colour = "#525e61",
        alpha = 0.8
    ) +
    scale_colour_manual(values = colour_values, drop = FALSE) +
    scale_fill_manual(values = fill_values, drop = FALSE) +
    scale_shape_manual(values = shape_values, drop = FALSE) +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y
    ) +
    theme_prism() +
    theme(
        plot.background = element_blank()
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output)
