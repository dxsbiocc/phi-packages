#!/usr/bin/env Rscript

# Template-ID: scatter-lollipop-polar
#
# Purpose:
#   Draw a polar lollipop chart of a numeric value around a full circle
#   of categories, with rotated labels.
#
# Inputs:
#   A table with one row per category. Default example:
#     - age_group: category label
#     - mean_prev: numeric value
#     - group: colour group
#
# Output:
#   A PDF, PNG, or SVG polar lollipop plot.
#
# Dependencies:
#   ggplot2, readr
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (column names and labels).
#   Edit DATA PREPARATION to change category order or label rotation.
#   Edit PLOT only when the chart geometry or styling must change.
#
# Scientific assumptions:
#   Each row is one category with a comparable numeric value.
#   Polar layout and label rotation are display choices.
#   Values are treated as already summarized.

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
        x = "age_group",
        y = "mean_prev",
        colour = "group"
    ),
    labels = list(
        title = "Polar Lollipop Plot",
        x = "",
        y = ""
    )
)

load_packages(c("ggplot2", "readr"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
require_columns(df, config$columns)

x_col <- config$columns$x
y_col <- config$columns$y
colour_col <- config$columns$colour

df[[x_col]] <- factor(df[[x_col]], levels = unique(df[[x_col]]))
df[[colour_col]] <- factor(df[[colour_col]], levels = unique(df[[colour_col]]))
df[[y_col]] <- as.numeric(df[[y_col]])

angle <- 90 - 360 * (seq_len(nrow(df)) - 0.5) / nrow(df)
df$label_angle <- ifelse(angle < -90, angle + 180, angle)
df$label_hjust <- ifelse(angle < -90, 1.2, -0.2)
df$value_angle <- ifelse(angle < -90, angle + 180, angle - 6)
df$value_hjust <- ifelse(angle < -90, 1.4, -0.4)
df$value_label <- round(df[[y_col]], 2)

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
colour_values <- c("#1aafd0", "#2dde98", "#ffc168")

p <- ggplot(df, aes(
    x = .data[[config$columns$x]],
    y = .data[[config$columns$y]],
    colour = .data[[config$columns$colour]]
)) +
    geom_hline(yintercept = c(0, 0.1, 0.2), colour = "grey90") +
    geom_segment(aes(
        x = .data[[config$columns$x]],
        xend = .data[[config$columns$x]],
        y = 0,
        yend = .data[[config$columns$y]]
    )) +
    geom_point(size = 5, shape = 21, fill = "#ffffff") +
    geom_point(size = 3) +
    geom_text(aes(
        label = .data[[config$columns$x]],
        angle = label_angle,
        hjust = label_hjust
    )) +
    geom_text(aes(
        label = value_label,
        angle = value_angle,
        hjust = value_hjust
    ), vjust = 2) +
    annotate(
        "rect",
        xmin = -Inf,
        xmax = Inf,
        ymin = -0.1,
        ymax = 0,
        fill = "#81947a",
        alpha = 0.5
    ) +
    annotate("text", x = 0, y = -0.1, label = "Age", size = 6, vjust = 0.5) +
    coord_polar() +
    scale_colour_manual(values = colour_values) +
    scale_y_continuous(limits = c(-0.1, 0.3)) +
    guides(colour = "none") +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y
    ) +
    theme(
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.background = element_blank(),
        plot.title = element_text(size = 14, hjust = 0.5)
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output)
