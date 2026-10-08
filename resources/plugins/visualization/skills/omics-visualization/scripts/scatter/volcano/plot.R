#!/usr/bin/env Rscript

# Template-ID: scatter-volcano
#
# Purpose:
#   Draw a volcano plot of log2 fold-change versus -log10(q-value), with
#   labels for the strongest up- and down-regulated genes.
#
# Inputs:
#   A table with one row per gene. Default example:
#     - log2FC: log2 fold-change
#     - qvalue: adjusted p-value
#     - group: direction label (None / Up / Down)
#     - symbol: gene label
#
# Output:
#   A PDF, PNG, or SVG volcano plot.
#
# Dependencies:
#   ggplot2, readr, ggprism, ggrepel
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (column names and labels).
#   Edit DATA PREPARATION to change the q-value floor, significance
#   cutoff, or how many labels are drawn per direction.
#   Edit PLOT only when the chart geometry or styling must change.
#
# Scientific assumptions:
#   Each row is one gene with a log2 fold-change and an adjusted p-value.
#   q-values are floored at 1e-30 so -log10(q) stays finite.
#   Labels are taken from genes with q-value < 0.05.

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
        x = "log2FC",
        y = "qvalue",
        colour = "group",
        label = "symbol"
    ),
    labels = list(
        title = "Yap1-KO",
        x = "log2FC",
        y = "-log10(p.adjust)"
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
df[[y_col]] <- pmax(as.numeric(df[[y_col]]), 1e-30)
df[[colour_col]] <- factor(df[[colour_col]], levels = unique(df[[colour_col]]))
df$neg_log10 <- -log10(df[[y_col]])

sig <- df[df[[y_col]] < 0.05, , drop = FALSE]
text_data <- do.call(rbind, lapply(
    split(sig, sig[[x_col]] < 0),
    function(sub) {
        sub[head(order(abs(sub[[x_col]]), decreasing = TRUE), 10), , drop = FALSE]
    }
))

grad_cols <- grDevices::colorRampPalette(c("#DFE1B8", "#FFFFFF", "#F3D7CE"))(1000)
grad <- matrix(grad_cols, nrow = 1)

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
colour_values <- c(
    Down = "#a3ad62",
    None = "grey70",
    Up = "#df91a3"
)

p <- ggplot(df, aes(
    x = .data[[config$columns$x]],
    y = neg_log10,
    colour = .data[[config$columns$colour]],
    label = .data[[config$columns$label]]
)) +
    annotation_raster(
        grad,
        xmin = -Inf,
        xmax = Inf,
        ymin = -Inf,
        ymax = Inf,
        interpolate = FALSE
    ) +
    geom_point(aes(size = abs(.data[[config$columns$x]])), alpha = 0.8, stroke = 0) +
    geom_text_repel(
        data = text_data,
        show.legend = FALSE,
        colour = "#000000"
    ) +
    geom_vline(xintercept = 0, linewidth = 1) +
    annotate(
        "segment",
        x = -0.2,
        y = c(10, 20, 30),
        xend = 0,
        yend = c(10, 20, 30),
        linewidth = 1
    ) +
    annotate(
        "text",
        x = -0.3,
        y = c(10, 20, 30),
        label = c(10, 20, 30),
        hjust = 1,
        colour = "black"
    ) +
    scale_colour_manual(values = colour_values, drop = FALSE) +
    scale_size_continuous(range = c(1, 8)) +
    scale_x_continuous(limits = c(-12, 12)) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
    guides(
        size = "none",
        colour = guide_legend(override.aes = list(size = 6))
    ) +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y
    ) +
    theme_prism() +
    theme(
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        plot.background = element_blank()
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output)
