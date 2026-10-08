#!/usr/bin/env Rscript

# Template-ID: heatmap-corr-dot
#
# Purpose:
#   Draw a rectangular association grid: colour encodes a supplied
#   signed coefficient, shape encodes a supplied p-value cutoff.
#
# Inputs:
#   One row per cell. Default example:
#     - cancer: column category (for example tumour type)
#     - cell: row category (for example immune subset)
#     - rho: signed association (Spearman, partial r, …)
#     - p: supplied p-value
#
# Output:
#   A PDF, PNG, or SVG correlation-dot heatmap.
#
# Dependencies:
#   ggplot2, readr, ggprism
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (columns, p_cutoff, labels).
#   Edit DATA PREPARATION to change axis order.
#   Edit PLOT only when the chart geometry or styling must change.
#
# Scientific assumptions:
#   rho and p are supplied results. This script does not correlate
#   expression with cell fractions, and does not adjust p-values.
#   Shape is a display of p versus p_cutoff, not a new test.
#   File order of cancer and cell is the display order (first y
#   level at the bottom, matching ggplot discrete defaults).

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
        x = "cancer",
        y = "cell",
        value = "rho",
        p = "p"
    ),
    p_cutoff = 0.05,
    point_size = 6,
    palette = "Diverging.TealRose",
    labels = list(
        title = "",
        x = NULL,
        y = NULL,
        colour = "Partial_Cor",
        shape = NULL
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
v_col <- config$columns$value
p_col <- config$columns$p
p_cutoff <- config$p_cutoff

df[[x_col]] <- factor(df[[x_col]], levels = unique(df[[x_col]]))
df[[y_col]] <- factor(df[[y_col]], levels = unique(df[[y_col]]))
df[[v_col]] <- as.numeric(df[[v_col]])
df[[p_col]] <- as.numeric(df[[p_col]])

ok <- is.finite(df[[v_col]]) & is.finite(df[[p_col]]) & df[[p_col]] >= 0
if (any(!ok)) {
    message(
        "Dropped ", sum(!ok),
        " row(s) with non-finite rho or invalid p."
    )
    df <- df[ok, , drop = FALSE]
}
if (!nrow(df)) {
    stop("No rows left to plot.", call. = FALSE)
}

sig_hi <- paste0("p \u2264 ", p_cutoff)
sig_lo <- paste0("p > ", p_cutoff)
df$sig <- factor(
    ifelse(df[[p_col]] <= p_cutoff, sig_hi, sig_lo),
    levels = c(sig_hi, sig_lo)
)

rho_max <- max(abs(df[[v_col]]), na.rm = TRUE)
if (!is.finite(rho_max) || rho_max == 0) {
    rho_max <- 1
}
ramp <- palette_colors(config$palette)
div_cols <- grDevices::colorRampPalette(
    c(ramp[[1]], "#ffffff", ramp[[length(ramp)]])
)(50)
shape_values <- c(15, 7)
names(shape_values) <- levels(df$sig)

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p <- ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    geom_point(
        aes(colour = .data[[v_col]], shape = sig),
        size = config$point_size,
        na.rm = TRUE
    ) +
    scale_shape_manual(
        name = config$labels$shape,
        values = shape_values,
        breaks = levels(df$sig),
        drop = FALSE
    ) +
    scale_colour_gradientn(
        name = config$labels$colour,
        colours = div_cols,
        limits = c(-rho_max, rho_max),
        breaks = round(c(-rho_max + 0.1, 0, rho_max - 0.1), 1)
    ) +
    scale_y_discrete(position = "right") +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y
    ) +
    theme_prism() +
    theme(
        panel.background = element_blank(),
        panel.border = element_rect(colour = "grey50", fill = NA, linewidth = 1),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text = element_text(size = 10),
        axis.ticks.length = unit(0, "cm"),
        axis.line = element_blank(),
        panel.grid = element_blank(),
        plot.background = element_blank(),
        legend.background = element_blank()
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output, width = 10, height = 5)
