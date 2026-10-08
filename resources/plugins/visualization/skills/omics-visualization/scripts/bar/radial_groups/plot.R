#!/usr/bin/env Rscript

# Template-ID: bar-radial-groups
#
# Purpose:
#   Draw a grouped radial bar chart with gene labels around the circle
#   and group arcs in the inner ring.
#
# Inputs:
#   A table with one row per gene-by-group. Default example:
#     - gene: category label
#     - group: group name
#     - count: numeric value
#
# Output:
#   A PDF, PNG, or SVG grouped radial bar chart.
#
# Dependencies:
#   ggplot2, readr, geomtextpath
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (column names and labels).
#   Edit DATA PREPARATION to change spacer rows, sort order, or label
#   rotation.
#   Edit PLOT only when the chart geometry or styling must change.
#
# Scientific assumptions:
#   Each row is one gene-by-group pair with a comparable numeric value.
#   Spacer rows between groups are visual, not additional observations.
#   Polar coordinates are a display choice, not a statistical transform.

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
        x = "gene",
        y = "count",
        group = "group"
    ),
    labels = list(
        title = "",
        x = "",
        y = ""
    )
)

load_packages(c("ggplot2", "readr", "geomtextpath"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
require_columns(df, config$columns)

x_col <- config$columns$x
y_col <- config$columns$y
group_col <- config$columns$group

df[[x_col]] <- as.character(df[[x_col]])
df[[y_col]] <- as.numeric(df[[y_col]])
df[[group_col]] <- as.character(df[[group_col]])
groups <- unique(df[[group_col]])

# One spacer row per group so arcs do not wrap into the next group.
spacer <- df[rep(NA_integer_, length(groups)), , drop = FALSE]
spacer[[x_col]] <- paste0(x_col, seq_along(groups))
spacer[[group_col]] <- groups
spacer[[y_col]] <- NA_real_
df <- rbind(df, spacer)

df[[group_col]] <- factor(df[[group_col]], levels = groups)
df <- df[order(
    df[[group_col]], is.na(df[[y_col]]), df[[x_col]]
), , drop = FALSE]
df$id <- seq_len(nrow(df))
df$angle <- 90 - 360 * (df$id - 0.5) / nrow(df)
df$text_angle <- ifelse(df$angle < -90, df$angle + 180, df$angle)
df$text_hjust <- ifelse(df$angle > -90, 0.3, 0.7)

base_anno <- do.call(rbind, lapply(groups, function(g) {
    ids <- df$id[df[[group_col]] == g]
    data.frame(
        group = g,
        start = min(ids),
        end = max(ids) - 1,
        mid = (min(ids) + max(ids) - 1) / 2,
        stringsAsFactors = FALSE
    )
}))

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
fill_values <- c(
    "#009392", "#39b185", "#9ccb86", "#e9e29c", "#eeb479", "#e88471", "#cf597e"
)

p <- ggplot(df, aes(
    x = id,
    y = .data[[y_col]],
    fill = .data[[group_col]]
)) +
    geom_col(position = position_dodge2(), show.legend = FALSE) +
    geom_text(
        aes(
            y = .data[[y_col]] + 18,
            label = .data[[x_col]],
            angle = text_angle,
            hjust = text_hjust
        ),
        size = 4,
        alpha = 0.6
    ) +
    geom_segment(
        data = base_anno,
        aes(x = start, xend = end, y = -5, yend = -5),
        color = "grey40",
        inherit.aes = FALSE
    ) +
    geom_textpath(
        data = base_anno,
        aes(x = mid, y = -18, label = group),
        color = "grey40",
        inherit.aes = FALSE,
        size = 6
    ) +
    scale_fill_manual(values = fill_values) +
    scale_y_continuous(limits = c(-100, 120)) +
    coord_polar() +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y
    ) +
    theme_minimal() +
    theme(
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        panel.grid = element_blank(),
        plot.background = element_blank()
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output)
