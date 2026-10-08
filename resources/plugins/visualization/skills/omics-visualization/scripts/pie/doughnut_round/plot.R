#!/usr/bin/env Rscript

# Template-ID: pie-doughnut-round
#
# Purpose:
#   Draw a doughnut (ring) chart whose slices are rounded rectangles,
#   matching the ECharts pie-borderRadius geometry.
#
# Inputs:
#   A two-column table (tsv/csv). Default example:
#     - name: slice label
#     - value: numeric size
#
# Output:
#   A PDF, PNG, or SVG rounded doughnut chart.
#
# Dependencies:
#   ggplot2, readr, ggprism, gground
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (column names and labels).
#   Edit DATA PREPARATION to change category order or reshape the table.
#   Edit PLOT only when the chart geometry or styling must change.
#
# Scientific assumptions:
#   Each row is one mutually exclusive part of a whole.
#   Inner hole, slice padding, and corner radius are display offsets,
#   not missing data.
#   Values are treated as already summarized, not as raw replicates.

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
        x = "name",
        y = "value"
    ),
    labels = list(
        title = "Rounded Doughnut Chart",
        x = "",
        y = ""
    )
)

load_packages(c("ggplot2", "readr", "ggprism", "gground"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
require_columns(df, config$columns)
x_col <- config$columns$x
y_col <- config$columns$y
df[[x_col]] <- factor(df[[x_col]], levels = unique(df[[x_col]]))
df[[y_col]] <- as.numeric(df[[y_col]])
df <- df[order(df[[x_col]]), , drop = FALSE]

# ECharts pie-borderRadius: radius ['40%', '70%'], padAngle 5, borderRadius 10.
inner_radius <- 0.40
outer_radius <- 0.70
pad_frac <- 5 / 360
corner_radius <- 10

n_slice <- nrow(df)
total <- sum(df[[y_col]])
available <- 1 - n_slice * pad_frac
if (available <= 0) {
    stop("padAngle leaves no angle for slices; reduce the gap.", call. = FALSE)
}
share <- df[[y_col]] / total * available
df$ymax <- cumsum(share + pad_frac)
df$ymin <- df$ymax - share
df$xmin <- inner_radius
df$xmax <- outer_radius

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p <- ggplot(df, aes(
    xmin = xmin,
    xmax = xmax,
    ymin = ymin,
    ymax = ymax,
    fill = .data[[x_col]]
)) +
    geom_round_rect(radius = corner_radius, colour = NA) +
    scale_fill_manual(
        values = palette_colors("Qualitative.Safe", n = nlevels(df[[x_col]]))
    ) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    coord_polar(theta = "y", start = 0, direction = 1, clip = "off") +
    labs(title = config$labels$title, x = NULL, y = NULL, fill = "") +
    theme_prism() +
    theme(
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        plot.background = element_blank(),
        legend.background = element_blank()
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output, width = 8, height = 8)
