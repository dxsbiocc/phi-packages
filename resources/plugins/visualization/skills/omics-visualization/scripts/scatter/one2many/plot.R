#!/usr/bin/env Rscript

# Template-ID: scatter-one2many
#
# Purpose:
#   Draw a one-to-many scatter plot of correlations around a central
#   feature, with Bezier links to each related item.
#
# Inputs:
#   A table with one row per related item. Default example:
#     - label: item name
#     - category: grouping label
#     - correlation: signed correlation
#
# Output:
#   A PDF, PNG, or SVG one-to-many scatter plot.
#
# Dependencies:
#   ggplot2, readr, ggforce, ggnewscale
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (column names and labels).
#   Edit DATA PREPARATION to change how points are placed or how
#   Bezier control points are built.
#   Edit PLOT only when the chart geometry or styling must change.
#
# Scientific assumptions:
#   Each row is one related item with a signed correlation.
#   Left/right placement encodes correlation sign, not a spatial axis.
#   Bezier curves are a display link to the central feature.

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
        label = "label",
        group = "category",
        correlation = "correlation"
    ),
    labels = list(
        title = "One-to-Many Scatter Plot",
        x = "",
        y = ""
    )
)

load_packages(c("ggplot2", "readr", "ggforce", "ggnewscale"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
require_columns(df, config$columns)

label_col <- config$columns$label
group_col <- config$columns$group
corr_col <- config$columns$correlation

df[[corr_col]] <- as.numeric(df[[corr_col]])
df[[group_col]] <- factor(df[[group_col]], levels = unique(df[[group_col]]))
df <- df[order(df[[group_col]], df[[label_col]]), ]
df$x <- ifelse(df[[corr_col]] < 0, -1, 1)
df$y <- rev(seq_len(nrow(df)))
df$y <- df$y - mean(df$y)
df$abs_corr <- abs(df[[corr_col]])
df$text_hjust <- ifelse(df$x > 0, 0, 1)
df$text_nudge <- ifelse(df$x > 0, 0.2, -0.2)

adjust <- 0.5
bezier_data <- do.call(rbind, lapply(seq_len(nrow(df)), function(i) {
    if (df$x[[i]] < 0) {
        data.frame(
            x = c(0, -adjust, df$x[[i]] + adjust, df$x[[i]]),
            y = c(0, 0, df$y[[i]], df$y[[i]]),
            type = "neg",
            group = i
        )
    } else {
        data.frame(
            x = c(0, adjust, df$x[[i]] - adjust, df$x[[i]]),
            y = c(0, 0, df$y[[i]], df$y[[i]]),
            type = "pos",
            group = i
        )
    }
}))

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
link_colours <- c(neg = "#84afa2", pos = "#f29556")
fill_values <- c(
    "#d7689d", "#f3d337", "#a0c443", "#66b5ae", "#5da4dc", "#40769e"
)

p <- ggplot(df, aes(
    x = x,
    y = y,
    label = .data[[config$columns$label]],
    fill = .data[[config$columns$group]],
    size = abs_corr
)) +
    geom_bezier(
        data = bezier_data,
        aes(x = x, y = y, group = group, colour = type),
        linewidth = 0.4,
        inherit.aes = FALSE,
        show.legend = FALSE
    ) +
    scale_colour_manual(values = link_colours) +
    new_scale_colour() +
    annotate("point", x = 0, y = 0, size = 8, colour = "#2B393B") +
    annotate(
        "text",
        x = 1.5,
        y = 0,
        label = "ethyl glucuronide",
        size = 6,
        hjust = 0.5,
        vjust = 0.5,
        colour = "#2B393B",
        fontface = "italic"
    ) +
    geom_point(colour = "white", show.legend = FALSE) +
    geom_point(shape = 21) +
    geom_text(
        aes(hjust = text_hjust),
        nudge_x = df$text_nudge,
        vjust = 0.5,
        size = 4,
        alpha = 0.8,
        show.legend = FALSE
    ) +
    annotate(
        "segment",
        x = -0.5,
        y = -13,
        xend = -2.5,
        yend = -13,
        arrow = arrow(angle = 30, length = unit(0.03, "npc"), type = "closed")
    ) +
    annotate("text",
        x = -2.3, y = -14,
        label = "Negative correlation",
        hjust = 0,
        colour = "#2B393B"
    ) +
    annotate(
        "segment",
        x = 0.5,
        y = -13,
        xend = 2.5,
        yend = -13,
        arrow = arrow(angle = 20, length = unit(0.03, "npc"), type = "closed")
    ) +
    annotate("text",
        x = 2.3, y = -14,
        label = "Positive correlation",
        hjust = 1,
        colour = "#2B393B"
    ) +
    scale_fill_manual(
        name = "Category",
        values = fill_values,
        guide = guide_legend(override.aes = list(size = 3))
    ) +
    scale_size_continuous(name = "Correlation", range = c(1, 6)) +
    scale_x_continuous(limits = c(-3, 3)) +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y
    ) +
    theme(
        panel.background = element_blank(),
        plot.background = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.background = element_blank(),
        legend.position = "top"
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output)
