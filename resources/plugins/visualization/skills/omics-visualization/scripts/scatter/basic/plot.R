#!/usr/bin/env Rscript

# Template-ID: scatter-basic
#
# Purpose:
#   Draw a basic scatter plot of two numeric variables.
#   style = "rug" adds coloured marginal rugs on the bottom (x) and
#   left (y) axes.
#
# Inputs:
#   A two-column table (tsv/csv). Default example:
#     - x: first numeric axis
#     - y: second numeric axis
#   Rug example (example_rug.tsv): wt, mpg (mtcars).
#
# Output:
#   A PDF, PNG, or SVG scatter plot.
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
#   Grouped histogram rugs are scatter-hist-rug. Marginal boxplots
#   are scatter-marginal.
#
# Scientific assumptions:
#   Each row is one paired observation of two numeric variables.
#   Rugs are a display of the marginal sample; no test is run.

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
    # "basic": plain points (default preview).
    # "rug": points plus x/y marginal rugs in catalog colours.
    style = "basic",
    columns = list(
        x = "x",
        y = "y"
    ),
    labels = list(
        title = "Scatter Plot",
        x = "x",
        y = "y"
    ),
    rug = list(
        point = "Qualitative.Safe:12",
        x = "Qualitative.Bold:3",
        y = "Qualitative.Bold:7",
        point_size = 2.4,
        rug_linewidth = 1.1,
        rug_alpha = 0.85,
        rug_length = 0.03,
        base_size = 13,
        width = 6.2,
        height = 5
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
df[[x_col]] <- as.numeric(df[[x_col]])
df[[y_col]] <- as.numeric(df[[y_col]])

pal_ref <- function(ref) {
    parts <- strsplit(ref, ":", fixed = TRUE)[[1]]
    cols <- palette_colors(parts[[1]])
    if (length(parts) > 1L) cols[[as.integer(parts[[2]])]] else cols
}

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
if (identical(config$style, "rug")) {
    rc <- config$rug
    p <- ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]])) +
        geom_point(
            size = rc$point_size,
            colour = pal_ref(rc$point)
        ) +
        geom_rug(
            sides = "b",
            colour = pal_ref(rc$x),
            linewidth = rc$rug_linewidth,
            alpha = rc$rug_alpha,
            length = grid::unit(rc$rug_length, "npc")
        ) +
        geom_rug(
            sides = "l",
            colour = pal_ref(rc$y),
            linewidth = rc$rug_linewidth,
            alpha = rc$rug_alpha,
            length = grid::unit(rc$rug_length, "npc")
        ) +
        labs(
            title = config$labels$title,
            x = config$labels$x,
            y = config$labels$y
        ) +
        theme_bw(base_size = rc$base_size) +
        theme(
            panel.grid.minor = element_blank(),
            plot.background = element_rect(fill = "white", colour = NA),
            panel.background = element_rect(fill = "white", colour = NA),
            axis.text = element_text(colour = "black")
        )
    if (is.null(config$labels$title) ||
        !nzchar(as.character(config$labels$title))) {
        p <- p + theme(plot.title = element_blank())
    }
    save_ggplot(p, io$output, width = rc$width, height = rc$height)
    quit(save = "no")
}

p <- ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    geom_point(size = 2.5, colour = "#1aafd0") +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y
    ) +
    theme_prism() +
    theme(plot.background = element_blank())

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output)
