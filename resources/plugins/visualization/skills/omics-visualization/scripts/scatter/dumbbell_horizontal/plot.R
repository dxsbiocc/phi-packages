#!/usr/bin/env Rscript

# Template-ID: scatter-dumbbell-horizontal
#
# Purpose:
#   Draw a horizontal dumbbell plot comparing a numeric value across two
#   groups on a shared categorical axis.
#
# Inputs:
#   A table with one row per category-by-group. Default example:
#     - age_group: category label
#     - mean_prev: numeric value
#     - Sex: group (two levels)
#
# Output:
#   A PDF, PNG, or SVG horizontal dumbbell plot.
#
# Dependencies:
#   ggplot2, readr, ggprism, scales
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (column names and labels).
#   Edit DATA PREPARATION to change which group sits left or category
#   order.
#   Edit PLOT only when the chart geometry or styling must change.
#
# Scientific assumptions:
#   Each row is one category-by-group pair with a comparable numeric value.
#   Values are treated as already summarized, not as raw replicates.
#   Percent labels are a display transform of the same values.

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
        x = "mean_prev",
        y = "age_group",
        group = "Sex"
    ),
    labels = list(
        title = "Mean Prevalence of obesity",
        x = "",
        y = ""
    )
)

load_packages(c("ggplot2", "readr", "ggprism", "scales"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
require_columns(df, config$columns)

x_col <- config$columns$x
y_col <- config$columns$y
group_col <- config$columns$group

df[[y_col]] <- factor(df[[y_col]], levels = unique(df[[y_col]]))
df[[group_col]] <- factor(df[[group_col]], levels = unique(df[[group_col]]))
df[[x_col]] <- as.numeric(df[[x_col]])

groups <- levels(df[[group_col]])
df$label_hjust <- ifelse(df[[group_col]] == groups[[1]], 1.2, -0.2)
df$pct_label <- scales::percent(df[[x_col]], accuracy = 0.1)

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
colour_values <- c("#f784b6", "#79ceb8")

p <- ggplot(df, aes(
    x = .data[[config$columns$x]],
    y = .data[[config$columns$y]],
    colour = .data[[config$columns$group]]
)) +
    geom_line(aes(group = .data[[config$columns$y]]), colour = "#cba7a2") +
    geom_point(size = 5, shape = 21, fill = "#ffffff") +
    geom_point(size = 3) +
    geom_text(aes(label = pct_label, hjust = label_hjust)) +
    annotate("segment", x = 0.08, y = 16, xend = 0.27, yend = 16, linewidth = 1.2) +
    annotate(
        "text",
        x = 0.08,
        y = 16.7,
        label = "Male",
        colour = "#79ceb8",
        hjust = 0,
        size = 5
    ) +
    annotate(
        "text",
        x = 0.27,
        y = 16.7,
        label = "Female",
        colour = "#f784b6",
        hjust = 1,
        size = 5
    ) +
    scale_colour_manual(values = colour_values) +
    scale_x_continuous(limits = c(0.08, 0.27)) +
    scale_y_discrete(expand = expansion(mult = c(0.1, 0.15))) +
    guides(colour = "none", fill = "none") +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y
    ) +
    theme_prism() +
    theme(
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        plot.background = element_blank()
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output)
