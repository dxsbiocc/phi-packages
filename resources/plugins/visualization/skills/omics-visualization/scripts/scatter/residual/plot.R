#!/usr/bin/env Rscript

# Template-ID: scatter-residual
#
# Purpose:
#   Draw an x-y scatter with an OLS fit, residual stems from each
#   observation to its predicted value, open circles on the fit, and
#   observed points sized and filled by absolute residual.
#
# Inputs:
#   One row per observation. Default example (diamonds sample):
#     - carat, price
#   Optional columns predicted and/or abs_rd: if missing, they are
#   computed from lm(y ~ x) as a display residual (not a test).
#
# Output:
#   A PDF, PNG, or SVG residual scatter.
#
# Dependencies:
#   ggplot2, readr, ggprism
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   Edit CONFIG for columns, palette ends, theme, and sizes.
#   Grouped multi-fit equations stay scatter-fit-group. A plain
#   correlation label stays scatter-correlation. Do not run
#   multivariate models or formal residual diagnostics here.
#
# Scientific assumptions:
#   Each row is one paired observation. The line is OLS y ~ x.
#   Absolute residual is |y - predicted|; encoding size/fill is a
#   display choice, not an inference.

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
        x = "carat",
        y = "price",
        predicted = NULL,  # optional; else lm(y ~ x)
        abs_rd = NULL      # optional; else abs(y - predicted)
    ),
    labels = list(
        title = NULL,
        x = "carat",
        y = "price",
        fill = "abs_rd",
        size = "abs_rd"
    ),
    # Continuous fill from low residual to high residual.
    palettes = list(
        fill_low = "Qualitative.Paired:5",
        fill_high = "Qualitative.Paired:3"
    ),
    line_colour = "grey75",
    line_width = 1.0,
    segment_alpha = 0.2,
    segment_colour = "grey40",
    predicted_shape = 1,
    predicted_size = 1.6,
    predicted_colour = "grey20",
    point_shape = 21,
    point_colour = "black",
    point_stroke = 0.35,
    size_range = c(1.5, 8),
    theme = "prism",  # "prism" | "classic"
    size = list(width = 7.2, height = 5.4)
)

load_packages(c("ggplot2", "readr", "ggprism"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
pal_ref <- function(ref) {
    parts <- strsplit(ref, ":", fixed = TRUE)[[1]]
    cols <- palette_colors(parts[[1]])
    if (length(parts) > 1L) cols[[as.integer(parts[[2]])]] else cols
}

df <- as.data.frame(read_table_auto(io$input), stringsAsFactors = FALSE)
require_columns(df, list(config$columns$x, config$columns$y))
x_col <- config$columns$x
y_col <- config$columns$y
df[[x_col]] <- as.numeric(df[[x_col]])
df[[y_col]] <- as.numeric(df[[y_col]])
ok <- is.finite(df[[x_col]]) & is.finite(df[[y_col]])
if (any(!ok)) {
    message("Dropped ", sum(!ok), " non-finite row(s).")
    df <- df[ok, , drop = FALSE]
}
if (nrow(df) < 3L) {
    stop("Need at least 3 finite rows for a residual scatter.", call. = FALSE)
}

pred_col <- config$columns$predicted
abs_col <- config$columns$abs_rd
if (is.null(pred_col) || !pred_col %in% names(df)) {
    fit <- stats::lm(df[[y_col]] ~ df[[x_col]])
    df$predicted <- as.numeric(stats::predict(fit, newdata = df))
    pred_col <- "predicted"
} else {
    df[[pred_col]] <- as.numeric(df[[pred_col]])
}
if (is.null(abs_col) || !abs_col %in% names(df)) {
    df$abs_rd <- abs(df[[y_col]] - df[[pred_col]])
    abs_col <- "abs_rd"
} else {
    df[[abs_col]] <- as.numeric(df[[abs_col]])
}

fill_low <- pal_ref(config$palettes$fill_low)
fill_high <- pal_ref(config$palettes$fill_high)

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p <- ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    geom_segment(
        aes(xend = .data[[x_col]], yend = .data[[pred_col]]),
        alpha = config$segment_alpha,
        colour = config$segment_colour,
        linewidth = 0.3
    ) +
    geom_point(
        aes(fill = .data[[abs_col]], size = .data[[abs_col]]),
        shape = config$point_shape,
        colour = config$point_colour,
        stroke = config$point_stroke
    ) +
    geom_smooth(
        method = "lm",
        se = FALSE,
        colour = config$line_colour,
        linewidth = config$line_width,
        formula = y ~ x
    ) +
    geom_point(
        aes(y = .data[[pred_col]]),
        shape = config$predicted_shape,
        size = config$predicted_size,
        colour = config$predicted_colour
    ) +
    scale_fill_continuous(
        name = config$labels$fill,
        low = fill_low,
        high = fill_high
    ) +
    scale_size_continuous(
        name = config$labels$size,
        range = config$size_range
    ) +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y
    )

if (identical(config$theme, "prism")) {
    p <- p + theme_prism()
} else {
    p <- p + theme_classic()
}
p <- p + theme(
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    legend.background = element_rect(fill = "white", colour = NA)
)
if (is.null(config$labels$title) || !nzchar(as.character(config$labels$title))) {
    p <- p + theme(plot.title = element_blank())
}

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(
    p, io$output,
    width = config$size$width,
    height = config$size$height
)
