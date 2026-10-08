#!/usr/bin/env Rscript

# Template-ID: boxplot-polar
#
# Purpose:
#   Draw a radial violin plot of a numeric value across ordered
#   categories, with dashed radial grid lines and curved labels.
#
# Inputs:
#   A two-column table (tsv/csv). Default example:
#     - month: category label
#     - shap: numeric value
#
# Output:
#   A PDF, PNG, or SVG polar violin plot.
#
# Dependencies:
#   ggplot2, readr, geomtextpath
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (column names and labels).
#   Edit DATA PREPARATION to change category order or radial limits.
#   Edit PLOT only when the chart geometry or styling must change.
#
# Scientific assumptions:
#   Each row is one observation of a numeric value in a category.
#   Polar coordinates and grid segments are display transforms.
#   Category order in the file is the angular order.

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
        x = "month",
        y = "shap"
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
df[[x_col]] <- factor(df[[x_col]], levels = unique(df[[x_col]]))
df[[y_col]] <- as.numeric(df[[y_col]])

x_unique <- levels(df[[x_col]])
x_size <- length(x_unique)
ymin <- floor(min(df[[y_col]], na.rm = TRUE))
ymax <- ceiling(max(df[[y_col]], na.rm = TRUE))

grid_data <- data.frame(
    ymin = ymin,
    ymax = ymax
)
grid_data <- grid_data[rep(1, x_size), , drop = FALSE]
grid_data[[x_col]] <- x_unique

label_data <- data.frame(y = rep(ymax, x_size))
label_data[[x_col]] <- x_unique

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
fill_base <- c("#ed6ca4", "#fbb05b", "#acd372", "#7bc4e2")
fill_values <- if (length(fill_base) < x_size) {
    colorRampPalette(fill_base)(x_size)
} else {
    fill_base
}

p <- ggplot(df, aes(
    x = .data[[x_col]],
    y = .data[[y_col]],
    fill = .data[[x_col]]
)) +
    geom_segment(
        data = grid_data,
        aes(
            xend = .data[[x_col]],
            y = ymin,
            yend = ymax
        ),
        linetype = "dashed",
        colour = "grey50"
    ) +
    geom_violin(show.legend = FALSE) +
    annotate(
        "segment",
        x = c(0, x_size),
        xend = c(0, x_size),
        y = ymin,
        yend = ymax
    ) +
    geom_hline(yintercept = c(ymin, ymax)) +
    geom_textpath(
        data = label_data,
        aes(x = .data[[x_col]], y = y, label = .data[[x_col]]),
        vjust = 1.5,
        inherit.aes = FALSE
    ) +
    scale_fill_manual(values = fill_values) +
    scale_y_continuous(
        expand = c(0, 0),
        limits = c(-2, ymax + 1),
        breaks = seq(ymin, ymax, length.out = 5),
        guide = guide_axis(angle = -5)
    ) +
    coord_radial(start = -1.047198, end = 1.047198) +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y
    ) +
    theme(
        plot.background = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(face = "bold"),
        axis.ticks.x = element_blank(),
        axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank()
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output)
