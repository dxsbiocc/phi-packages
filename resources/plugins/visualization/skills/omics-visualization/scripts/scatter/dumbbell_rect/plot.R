#!/usr/bin/env Rscript

# Template-ID: scatter-dumbbell-rect
#
# Purpose:
#   Draw a rectangular dumbbell plot comparing two groups on a categorical
#   axis, with a rectangle spanning the two values in each category.
#
# Inputs:
#   A table with one row per category-by-group. Default example:
#     - age_group: category label
#     - mean_prev: numeric value
#     - Sex: group (two levels)
#
# Output:
#   A PDF, PNG, or SVG rectangular dumbbell plot.
#
# Dependencies:
#   ggplot2, readr, ggprism
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (column names and labels).
#   Edit DATA PREPARATION to change how the two group values are
#   reshaped into rectangles, or category order.
#   Edit PLOT only when the chart geometry or styling must change.
#
# Scientific assumptions:
#   Each row is one category-by-group pair with a comparable numeric value.
#   The rectangle joins the two group values in a category; it is not a
#   confidence interval.
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
        group = "Sex"
    ),
    labels = list(
        title = "",
        x = "",
        y = "Mean Prevalence of obesity"
    )
)

load_packages(c("ggplot2", "readr", "ggprism"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
require_columns(df, config$columns)

x_col <- config$columns$x
y_col <- config$columns$y
group_col <- config$columns$group

df[[x_col]] <- factor(df[[x_col]], levels = unique(df[[x_col]]))
df[[group_col]] <- factor(df[[group_col]], levels = unique(df[[group_col]]))
df[[y_col]] <- as.numeric(df[[y_col]])

group_levels <- levels(df[[group_col]])
rect_data <- data.frame(x = levels(df[[x_col]]), stringsAsFactors = FALSE)
names(rect_data)[1] <- x_col
for (g in group_levels) {
    sub <- df[df[[group_col]] == g, ]
    rect_data[[g]] <- sub[[y_col]][match(rect_data[[x_col]], as.character(sub[[x_col]]))]
}
rect_data[[x_col]] <- factor(rect_data[[x_col]], levels = levels(df[[x_col]]))
rect_data$xmin <- as.integer(rect_data[[x_col]]) + 0.1
rect_data$xmax <- as.integer(rect_data[[x_col]]) - 0.1
ymin_col <- group_levels[[1]]
ymax_col <- group_levels[[2]]

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
colour_values <- c("#f784b6", "#79ceb8")

p <- ggplot(df, aes(
    x = .data[[config$columns$x]],
    y = .data[[config$columns$y]],
    colour = .data[[config$columns$group]]
)) +
    geom_rect(
        data = rect_data,
        aes(
            xmin = xmin,
            xmax = xmax,
            ymin = .data[[ymin_col]],
            ymax = .data[[ymax_col]]
        ),
        fill = "#f3ead8",
        colour = "grey50",
        linewidth = 0.4,
        inherit.aes = FALSE
    ) +
    geom_point(size = 5, shape = 21, stroke = 1, fill = "#ffffff") +
    geom_point(size = 3) +
    scale_colour_manual(values = colour_values) +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y
    ) +
    theme_prism() +
    theme(
        legend.position = "inside",
        legend.position.inside = c(0.8, 0.8),
        axis.line = element_blank(),
        axis.text.x = element_text(angle = 90),
        panel.border = element_rect(fill = "transparent", colour = "black", linewidth = 1),
        plot.background = element_blank()
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output)
