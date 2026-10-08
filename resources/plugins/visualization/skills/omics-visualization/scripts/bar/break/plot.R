#!/usr/bin/env Rscript

# Template-ID: bar-break
#
# Purpose:
#   Draw a grouped bar chart with a broken y-axis so small and large series remain readable.
#
# Inputs:
#   A table with one row per category-by-group. Default example:
#     - day: category label
#     - value: numeric value
#     - category: group (series spanning different magnitudes)
#
# Output:
#   A PDF, PNG, or SVG grouped bar chart with a y-axis break.
#
# Dependencies:
#   ggplot2, readr, patchwork, ggprism
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (column names and labels).
#   Edit DATA PREPARATION to change the axis-break limits.
#   Edit PLOT only when the chart geometry or styling must change.
#
# Scientific assumptions:
#   Each row is one category-by-group pair with a comparable numeric value.
#   The y-axis break is a display device, not a statistical transform.
#   Panels are one zoomed window per magnitude cluster so within-group
#   differences stay visible.

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
        x = "day",
        y = "value",
        fill = "category"
    ),
    labels = list(
        title = "Bar Chart with Axis Break",
        x = "Day",
        y = "Value"
    )
)

load_packages(c("ggplot2", "readr", "patchwork", "ggprism"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
require_columns(df, config$columns)
x_col <- config$columns$x
y_col <- config$columns$y
fill_col <- config$columns$fill
df[[x_col]] <- factor(df[[x_col]], levels = unique(df[[x_col]]))
df[[fill_col]] <- factor(df[[fill_col]], levels = unique(df[[fill_col]]))
df[[y_col]] <- as.numeric(df[[y_col]])

# One zoomed y-window per magnitude cluster (high to low).
med <- tapply(df[[y_col]], df[[fill_col]], stats::median, na.rm = TRUE)
cluster_id <- setNames(as.integer(round(log10(pmax(med, 1)))), names(med))
axis_windows <- lapply(
    sort(unique(cluster_id), decreasing = TRUE),
    function(cl) {
        groups <- names(cluster_id)[cluster_id == cl]
        vals <- df[[y_col]][df[[fill_col]] %in% groups]
        ymin <- min(vals, na.rm = TRUE)
        ymax <- max(vals, na.rm = TRUE)
        span <- ymax - ymin
        if (!is.finite(span) || span <= 0) {
            span <- max(abs(ymax), 1)
        }
        pad <- span * 0.2
        lower <- if (ymin <= 0.25 * ymax) 0 else ymin - pad
        c(lower, ymax + pad)
    }
)

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
make_bars <- function(ylim, show_x) {
    ggplot(df, aes(
        x = .data[[x_col]],
        y = .data[[y_col]],
        fill = .data[[fill_col]]
    )) +
        geom_col(position = position_dodge(width = 0.85), width = 0.8) +
        scale_fill_manual(
            values = expand_palette(
                c("#d7689d", "#f3d337", "#a0c443", "#66b5ae"),
                nlevels(df[[fill_col]])
            )
        ) +
        coord_cartesian(ylim = ylim, expand = FALSE) +
        labs(x = if (show_x) config$labels$x else NULL, y = NULL, fill = "") +
        theme_prism() +
        theme(
            axis.text.x = if (show_x) element_text() else element_blank(),
            axis.ticks.x = if (show_x) element_line() else element_blank(),
            plot.background = element_blank()
        )
}

panels <- lapply(seq_along(axis_windows), function(i) {
    make_bars(axis_windows[[i]], i == length(axis_windows))
})
p <- wrap_plots(panels, ncol = 1, guides = "collect") +
    plot_annotation(title = config$labels$title)

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output)
