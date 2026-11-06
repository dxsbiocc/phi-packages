#!/usr/bin/env Rscript

# Template-ID: scatter-jitter
#
# Purpose:
#   Draw a jittered scatter plot of a numeric value across categories.
#   style = "dotplot" stacks observations into Wilkinson-style bins
#   along y (geom_dotplot), instead of random jitter.
#
# Inputs:
#   A two-column table (tsv/csv). Default example:
#     - day: category label
#     - value: numeric value
#   Dotplot example (example_dotplot.tsv):
#     - dose, len
#
# Output:
#   A PDF, PNG, or SVG jitter or stacked-dot plot.
#
# Dependencies:
#   ggplot2, readr, ggprism
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (column names and labels).
#   Edit DATA PREPARATION to change category order.
#   Edit PLOT only when the chart geometry or styling must change.
#   Enrichment term dotplots stay bar-enrichment-dot. Beeswarm
#   packing stays scatter-beeswarm / boxplot-basic style = beeswarm.
#
# Scientific assumptions:
#   Each row is one observation of a numeric value in a category.
#   Jitter and dot-binning are display transforms to reduce
#   overplotting, not statistical transforms of the values.

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
    # "jitter": random horizontal jitter (default preview).
    # "dotplot": Wilkinson stacked dots along y (geom_dotplot).
    style = "jitter",
    columns = list(
        x = "day",
        y = "value"
    ),
    labels = list(
        title = "Jitter Plot of Values by Day",
        x = "Value",
        y = "Day"
    ),
    dotplot = list(
        # Palette family expanded one colour per x level, or a single
        # "Family.Name:index" for a uniform fill.
        fill = "Qualitative.Bold",
        colour = "grey15",
        show_legend = FALSE,
        binwidth = 1.5,        # absolute y units; set NULL to let ggplot choose
        stackdir = "center",   # "center" | "up" | "down" | "centerwhole"
        stackratio = 1,
        dotsize = 1.15,
        stroke = 0.35,
        alpha = 1,
        base_size = 13,
        grid = TRUE,
        width = 6.5,
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
df[[x_col]] <- factor(df[[x_col]], levels = unique(df[[x_col]]))
df[[y_col]] <- as.numeric(df[[y_col]])

pal_ref <- function(ref) {
    parts <- strsplit(ref, ":", fixed = TRUE)[[1]]
    cols <- palette_colors(parts[[1]])
    if (length(parts) > 1L) cols[[as.integer(parts[[2]])]] else cols
}

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
if (identical(config$style, "dotplot")) {
    dc <- config$dotplot
    n_g <- nlevels(df[[x_col]])
    by_group <- !grepl(":", dc$fill, fixed = TRUE)
    if (by_group) {
        fill_vals <- expand_palette(palette_colors(dc$fill), n_g)
        names(fill_vals) <- levels(df[[x_col]])
        p <- ggplot(df, aes(
            x = .data[[x_col]],
            y = .data[[y_col]],
            fill = .data[[x_col]]
        ))
    } else {
        fill_vals <- pal_ref(dc$fill)
        p <- ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]]))
    }
    dp_args <- list(
        binaxis = "y",
        stackdir = dc$stackdir,
        stackratio = dc$stackratio,
        dotsize = dc$dotsize,
        stroke = dc$stroke,
        alpha = dc$alpha,
        colour = dc$colour
    )
    if (!by_group) {
        dp_args$fill <- fill_vals
    }
    if (!is.null(dc$binwidth)) {
        dp_args$binwidth <- dc$binwidth
    }
    p <- p + do.call(geom_dotplot, dp_args)
    if (by_group) {
        p <- p + scale_fill_manual(
            values = fill_vals,
            guide = if (isTRUE(dc$show_legend)) "legend" else "none"
        )
    }
    p <- p +
        scale_y_continuous(
            name = config$labels$y,
            expand = expansion(mult = c(0.02, 0.05))
        ) +
        labs(
            title = if (is.null(config$labels$title) ||
                !nzchar(as.character(config$labels$title))) {
                NULL
            } else {
                config$labels$title
            },
            x = config$labels$x
        ) +
        theme_classic(base_size = dc$base_size) +
        theme(
            axis.text = element_text(colour = "black"),
            axis.title = element_text(colour = "black"),
            axis.ticks = element_line(colour = "black"),
            panel.background = element_rect(fill = "white", colour = NA),
            plot.background = element_rect(fill = "white", colour = NA),
            panel.grid.major = if (isTRUE(dc$grid)) {
                element_line(colour = "grey88", linewidth = 0.3)
            } else {
                element_blank()
            },
            panel.grid.minor = element_blank()
        )
    if (is.null(config$labels$title) || !nzchar(as.character(config$labels$title))) {
        p <- p + theme(plot.title = element_blank())
    }
    save_ggplot(p, io$output, width = dc$width, height = dc$height)
    quit(save = "no")
}

colour_values <- c(
    "#8F499C", "#4185BE", "#6DC067", "#F6DB35",
    "#F78822", "#232121", "#CF202E"
)

p <- ggplot(df, aes(
    x = .data[[x_col]],
    y = .data[[y_col]]
)) +
    geom_jitter(
        aes(colour = .data[[x_col]]),
        size = 3,
        width = 0.18,
        alpha = 0.8,
        show.legend = FALSE
    ) +
    scale_colour_manual(values = colour_values) +
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
