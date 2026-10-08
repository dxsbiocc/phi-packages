#!/usr/bin/env Rscript

# Template-ID: scatter-beeswarm-group
#
# Purpose:
#   Draw a grouped beeswarm plot of log2 fold-change across contrasts,
#   with point size from -log10(adjusted p-value) and labels for
#   extreme genes.
#
# Inputs:
#   A table with one row per gene-by-contrast. Default example:
#     - cate: contrast label
#     - logFC: log2 fold-change
#     - sign: direction (Up / Down)
#     - adj.P.Val: adjusted p-value
#     - symbol: gene label
#
# Output:
#   A PDF, PNG, or SVG grouped beeswarm plot.
#
# Dependencies:
#   ggplot2, readr, ggprism, ggbeeswarm, ggrepel, ggnewscale
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (column names and labels).
#   Edit DATA PREPARATION to change contrast order or how many labels
#   are drawn per group.
#   Edit PLOT only when the chart geometry or styling must change.
#
# Scientific assumptions:
#   Each row is one gene in one contrast with a fold-change and p-value.
#   -log10(adj.P.Val) is a display transform of the adjusted p-value.
#   Beeswarm layout is a display choice to reduce overplotting.

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
        x = "cate",
        y = "logFC",
        group = "sign",
        pvalue = "adj.P.Val",
        label = "symbol"
    ),
    labels = list(
        title = "",
        x = "",
        y = "log2FC"
    )
)

load_packages(c("ggplot2", "readr", "ggprism", "ggbeeswarm", "ggrepel", "ggnewscale"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
require_columns(df, config$columns)

x_col <- config$columns$x
y_col <- config$columns$y
group_col <- config$columns$group
p_col <- config$columns$pvalue
label_col <- config$columns$label

keep <- c(x_col, y_col, group_col, p_col, label_col)
df <- df[stats::complete.cases(df[, keep]), ]
df[[x_col]] <- factor(df[[x_col]], levels = unique(df[[x_col]]))
df[[group_col]] <- factor(df[[group_col]], levels = unique(df[[group_col]]))
df[[y_col]] <- as.numeric(df[[y_col]])
df[[p_col]] <- pmax(as.numeric(df[[p_col]]), 1e-30)
df$neg_log10 <- -log10(df[[p_col]])

tile_data <- unique(df[, x_col, drop = FALSE])
text_data <- do.call(rbind, lapply(
    split(df, list(df[[x_col]], df[[group_col]]), drop = TRUE),
    function(sub) {
        sub[head(order(abs(sub[[y_col]]), decreasing = TRUE), 5), , drop = FALSE]
    }
))

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
sign_colour <- c(Down = "#2dde98", Up = "#ff6c5f")
cate_fill <- c("#a8dbd9", "#85c4c9", "#68abb8", "#4f90a6", "#3b738f")
cate_colour <- c("#595d7a", "#595d7a", "#595d7a", "#c5c5c5", "#c5c5c5")

p <- ggplot(df, aes(
    x = .data[[config$columns$x]],
    y = .data[[config$columns$y]],
    colour = .data[[config$columns$group]]
)) +
    geom_tile(
        data = tile_data,
        aes(
            x = .data[[config$columns$x]],
            y = 0,
            fill = .data[[config$columns$x]]
        ),
        height = 1.5,
        colour = "black",
        width = 0.9,
        show.legend = FALSE,
        inherit.aes = FALSE
    ) +
    geom_quasirandom(aes(size = neg_log10), width = 0.4, alpha = 0.5) +
    scale_colour_manual(
        name = "Change",
        values = sign_colour,
        labels = c("Down", "Up")
    ) +
    new_scale_colour() +
    geom_text(
        data = tile_data,
        aes(
            x = .data[[config$columns$x]],
            y = 0,
            label = .data[[config$columns$x]],
            colour = .data[[config$columns$x]]
        ),
        size = 4,
        show.legend = FALSE,
        inherit.aes = FALSE
    ) +
    geom_text_repel(
        data = text_data,
        aes(
            x = .data[[config$columns$x]],
            y = .data[[config$columns$y]],
            label = .data[[config$columns$label]]
        ),
        force = 2,
        size = 3,
        arrow = arrow(length = unit(0.008, "npc"), type = "open", ends = "last")
    ) +
    scale_fill_manual(values = cate_fill) +
    scale_colour_manual(values = cate_colour) +
    scale_size_binned(range = c(1, 4)) +
    guides(size = guide_bins(show.limits = TRUE)) +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y
    ) +
    theme_prism() +
    theme(
        plot.background = element_blank(),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 13, colour = "black", face = "bold"),
        axis.line.y = element_line(colour = "black", linewidth = 1),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank()
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output)
