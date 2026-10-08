#!/usr/bin/env Rscript

# Template-ID: heatmap-corr-triangle
#
# Purpose:
#   Draw a rectangular association grid whose cells are split on
#   the diagonal: one triangle is a supplied coefficient, the other
#   is a supplied p-value. Optional dots mark p-value cuts.
#
# Inputs:
#   One row per cell. Default example:
#     - cell: column category
#     - gene: row category
#     - cor: supplied coefficient
#     - p: supplied p-value
#
# Output:
#   A PDF, PNG, or SVG split-triangle heatmap.
#
# Dependencies:
#   ggplot2, readr, ggprism, ggnewscale
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (columns, which triangle is
#   |rho| vs p, p-dots, labels). Which triangle gets which column,
#   abs(rho), and dots on or off stay on this id. That is not
#   heatmap-corr-dot (shape) or heatmap-corr-bubble (size and
#   stars). Edit DATA PREPARATION to change axis order. Edit PLOT
#   only when the glyph must change.
#
# Scientific assumptions:
#   cor and p are supplied results. This script does not compute
#   Spearman, adjust p-values, or impute missing pairs. Triangle
#   vertices are a display of the categorical grid. Dots are a
#   display of supplied p versus CONFIG cuts.

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
        x = "cell",
        y = "gene",
        upper = "cor",
        lower = "p"
    ),
    upper_abs = TRUE,
    show_p_dots = TRUE,
    p_dots = c(
        "0.05" = 0.05,
        "0.01" = 0.01,
        "0.001" = 0.001
    ),
    labels = list(
        title = "",
        x = NULL,
        y = NULL,
        upper = "Upper triangle, |\u03c1|",
        lower = "Lower triangle, p"
    )
)

load_packages(c("ggplot2", "readr", "ggprism", "ggnewscale"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
need <- config$columns[c("x", "y", "upper", "lower")]
require_columns(df, need)

x_col <- config$columns$x
y_col <- config$columns$y
u_col <- config$columns$upper
l_col <- config$columns$lower

x_lv <- unique(as.character(df[[x_col]]))
y_lv <- unique(as.character(df[[y_col]]))
df[[x_col]] <- factor(df[[x_col]], levels = x_lv)
# File order of y is top-to-bottom.
df[[y_col]] <- factor(df[[y_col]], levels = rev(y_lv))
df[[u_col]] <- as.numeric(df[[u_col]])
df[[l_col]] <- as.numeric(df[[l_col]])

ok <- is.finite(df[[u_col]]) & is.finite(df[[l_col]]) & df[[l_col]] >= 0
if (any(!ok)) {
    message(
        "Dropped ", sum(!ok),
        " row(s) with non-finite coefficient or invalid p."
    )
    df <- df[ok, , drop = FALSE]
}
if (!nrow(df)) {
    stop("No rows left to plot.", call. = FALSE)
}

df$xi <- as.integer(df[[x_col]])
df$yi <- as.integer(df[[y_col]])
df$upper_fill <- df[[u_col]]
if (isTRUE(config$upper_abs)) {
    df$upper_fill <- abs(df$upper_fill)
}

n <- nrow(df)
gid <- as.character(seq_len(n))
upper <- data.frame(
    x = rep(df$xi, each = 3L) + c(0, 0, 1),
    y = rep(df$yi, each = 3L) + c(0, 1, 1),
    group = rep(gid, each = 3L),
    fill = rep(df$upper_fill, each = 3L)
)
lower <- data.frame(
    x = rep(df$xi, each = 3L) + c(0, 1, 1),
    y = rep(df$yi, each = 3L) + c(0, 0, 1),
    group = rep(gid, each = 3L),
    fill = rep(df[[l_col]], each = 3L)
)

pts <- NULL
if (isTRUE(config$show_p_dots)) {
    cuts <- unlist(config$p_dots)
    yoff <- seq(0.12, 0.52, length.out = max(length(cuts), 1L))
    rows <- lapply(seq_len(n), function(i) {
        hit <- df[[l_col]][[i]] < cuts
        if (!any(hit)) {
            return(NULL)
        }
        data.frame(
            x = df$xi[[i]] + 0.88,
            y = df$yi[[i]] + yoff[hit]
        )
    })
    pts <- do.call(rbind, rows)
}

u_max <- max(df$upper_fill, na.rm = TRUE)
if (!is.finite(u_max) || u_max <= 0) {
    u_max <- 1
}
p_max <- max(df[[l_col]], na.rm = TRUE)
if (!is.finite(p_max) || p_max <= 0) {
    p_max <- 1
}

# Project recommended pair, not yellow–blue: Sunset (gold–rose–purple)
# for |ρ|, BluGrn (mint–teal) for p. Low p stays light so black dots read.
upper_cols <- palette_colors("Quantitative.Sunset")
lower_cols <- palette_colors("Quantitative.BluGrn")
outline <- palette_colors("Qualitative.Safe")[[12]]
dot_col <- palette_colors("Brand.adidas")[[1]]

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p <- ggplot() +
    geom_polygon(
        data = upper,
        aes(x = x, y = y, group = group, fill = fill),
        colour = outline,
        linewidth = 0.25
    ) +
    scale_fill_gradientn(
        name = config$labels$upper,
        colours = upper_cols,
        limits = c(0, u_max),
        na.value = NA,
        guide = guide_colorbar(
            title.position = "top",
            title.hjust = 0,
            order = 1,
            barheight = unit(2.4, "cm")
        )
    ) +
    ggnewscale::new_scale_fill() +
    geom_polygon(
        data = lower,
        aes(x = x, y = y, group = group, fill = fill),
        colour = outline,
        linewidth = 0.25
    ) +
    scale_fill_gradientn(
        name = config$labels$lower,
        colours = lower_cols,
        limits = c(0, min(p_max, 1)),
        na.value = NA,
        guide = guide_colorbar(
            title.position = "top",
            title.hjust = 0,
            order = 2,
            barheight = unit(2.4, "cm")
        )
    ) +
    scale_x_continuous(
        breaks = seq_along(x_lv) + 0.5,
        labels = x_lv,
        expand = c(0, 0),
        limits = c(1, length(x_lv) + 1)
    ) +
    scale_y_continuous(
        breaks = seq_along(levels(df[[y_col]])) + 0.5,
        labels = levels(df[[y_col]]),
        expand = c(0, 0),
        limits = c(1, length(y_lv) + 1)
    ) +
    coord_equal(clip = "off") +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y
    ) +
    theme_prism() +
    theme(
        panel.grid = element_blank(),
        plot.background = element_blank(),
        panel.background = element_blank(),
        legend.background = element_blank(),
        panel.border = element_rect(colour = outline, fill = NA, linewidth = 0.55),
        axis.line = element_blank(),
        axis.ticks.length = unit(0, "pt"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.key = element_blank(),
        legend.title = element_text(size = 9, face = "bold"),
        legend.spacing.y = unit(8, "pt")
    )

if (!is.null(pts) && nrow(pts)) {
    p <- p + geom_point(
        data = pts,
        aes(x = x, y = y),
        size = 1.15,
        colour = dot_col,
        inherit.aes = FALSE
    )
}

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output, width = 9.6, height = 7.2)
