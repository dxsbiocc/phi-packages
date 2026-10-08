#!/usr/bin/env Rscript

# Template-ID: scatter-bubble
#
# Purpose:
#   Draw a bubble plot of GDP versus life expectancy, sized by population.
#
# Inputs:
#   A table with one row per country-year. Default example:
#     - GDP: x-axis (log-scaled for display)
#     - life: y-axis
#     - year: facet
#     - population: bubble size
#     - country: label
#
# Output:
#   A PDF, PNG, or SVG bubble chart.
#
# Dependencies:
#   ggplot2, readr, ggprism
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
#   Each row is one country at one year.
#   Log10 GDP is a display transform so small and large economies remain visible.
#   Bubble area encodes population; it is not a statistical weight.

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
        x = "GDP",
        y = "life",
        size = "population",
        colour = "country",
        facet = "year"
    ),
    labels = list(
        title = "GDP vs Life Expectancy",
        x = "GDP per capita",
        y = "Life expectancy"
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
size_col <- config$columns$size
colour_col <- config$columns$colour
facet_col <- config$columns$facet
df[[x_col]] <- as.numeric(df[[x_col]])
df[[y_col]] <- as.numeric(df[[y_col]])
df[[size_col]] <- as.numeric(df[[size_col]])
df[[facet_col]] <- factor(df[[facet_col]], levels = unique(df[[facet_col]]))

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p <- ggplot(df, aes(
    x = .data[[x_col]],
    y = .data[[y_col]],
    size = .data[[size_col]],
    colour = .data[[facet_col]]
)) +
    geom_point(alpha = 0.7) +
    scale_x_log10() +
    scale_size_area(max_size = 12) +
    scale_colour_manual(
        values = expand_palette(
            c("#2c98a0",
            "#b9257a"),
            nlevels(df[[facet_col]])
        )
    ) +
    facet_wrap(vars(.data[[facet_col]])) +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y,
        size = "Population",
        colour = "Year"
    ) +
    theme_prism() +
    theme(
        legend.position = "bottom",
        plot.background = element_blank()
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output, width = 10, height = 7)
