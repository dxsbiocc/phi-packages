#!/usr/bin/env Rscript

# Template-ID: heatmap-corr-bubble
#
# Purpose:
#   Draw a rectangular correlation bubble grid: fill encodes a supplied
#   signed coefficient, point size encodes magnitude (default -log10(p)
#   or a supplied size column), and stars mark supplied significance.
#
# Inputs:
#   One row per cell. Default example:
#     - CNV: column category
#     - Feature: row category
#     - rho: signed association (Spearman, partial r, …)
#     - p: supplied p-value
#     - Group, Type: optional facet columns
#
# Output:
#   A PDF, PNG, or SVG correlation-bubble heatmap.
#
# Dependencies:
#   ggplot2, readr, ggprism, scales, gtable
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (columns, fill limits, vline
#   intercepts, labels including x_size/y_size). Map size to a
#   supplied numeric column to encode abundance instead of -log10(p). Row-group bars come from
#   row_facet; extra column cuts use vline_intercepts. Facet vs
#   geom_vline is CONFIG, not a new id. Edit DATA PREPARATION to
#   change axis order. Edit PLOT only when the glyph must change.
#
# Scientific assumptions:
#   rho, p, and optional size are supplied results. This script does
#   not compute Spearman, adjust p-values, or estimate abundance.
#   Stars and -log10(p) are display transforms of the supplied p.

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
        x = "CNV",
        y = "Feature",
        value = "rho",
        p = "p",
        size = NULL,
        row_facet = "Type",
        col_facet = "Group"
    ),
    p_stars = c(
        "****" = 1e-4,
        "***" = 1e-3,
        "**" = 1e-2,
        "*" = 5e-2
    ),
    fill_limits = c(-0.55, 0.55),
    vline_intercepts = NULL,
    labels = list(
        title = "",
        x = NULL,
        y = NULL,
        fill = "Spearman's \u03c1",
        size = "-\u202flog10(P)",
        x_size = 20,
        y_size = 20
    )
)

load_packages(c("ggplot2", "readr", "ggprism", "scales", "gtable"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
need <- config$columns[c("x", "y", "value", "p")]
require_columns(df, need)

x_col <- config$columns$x
y_col <- config$columns$y
v_col <- config$columns$value
p_col <- config$columns$p
size_col <- config$columns$size
row_f <- config$columns$row_facet
col_f <- config$columns$col_facet

df[[x_col]] <- factor(df[[x_col]], levels = unique(as.character(df[[x_col]])))
# File order of y is top-to-bottom (ggplot's first level sits at the bottom).
y_levels <- unique(as.character(df[[y_col]]))
df[[y_col]] <- factor(df[[y_col]], levels = rev(y_levels))
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

if (!is.null(row_f) && row_f %in% names(df)) {
    df[[row_f]] <- factor(df[[row_f]], levels = unique(as.character(df[[row_f]])))
} else {
    row_f <- NULL
}
if (!is.null(col_f) && col_f %in% names(df)) {
    df[[col_f]] <- factor(df[[col_f]], levels = unique(as.character(df[[col_f]])))
} else {
    col_f <- NULL
}

if (!is.null(size_col) && size_col %in% names(df)) {
    df$size_enc <- as.numeric(df[[size_col]])
} else {
    df$size_enc <- -log10(pmax(df[[p_col]], .Machine$double.xmin))
    size_ok <- is.finite(df$size_enc)
    if (any(!size_ok)) {
        message("Dropped ", sum(!size_ok), " row(s) with non-finite size.")
        df <- df[size_ok, , drop = FALSE]
    }
}

star_cuts <- sort(unlist(config$p_stars), decreasing = TRUE)
star_lab <- function(p) {
    out <- character(length(p))
    for (i in seq_along(star_cuts)) {
        hit <- is.finite(p) & p < star_cuts[[i]]
        out[hit] <- names(star_cuts)[[i]]
    }
    out
}
df$star <- star_lab(df[[p_col]])

fill_lim <- as.numeric(config$fill_limits)
if (length(fill_lim) != 2L || any(!is.finite(fill_lim))) {
    abs_max <- max(abs(df[[v_col]]), na.rm = TRUE)
    if (!is.finite(abs_max) || abs_max == 0) {
        abs_max <- 1
    }
    fill_lim <- c(-abs_max, abs_max)
}

fill_cols <- palette_colors("Diverging.PRGn")
zero_at <- (0 - fill_lim[[1]]) / (fill_lim[[2]] - fill_lim[[1]])
zero_at <- min(max(zero_at, 0.05), 0.95)
n_fill <- length(fill_cols)
half <- (n_fill - 1L) / 2
fill_values <- c(
    seq(0, zero_at, length.out = half + 1L),
    seq(zero_at, 1, length.out = half + 1L)[-1]
)
outline <- palette_colors("Brand.Algolia")[[10]]

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p <- ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    geom_point(
        aes(fill = .data[[v_col]], size = size_enc),
        shape = 21,
        colour = outline,
        stroke = 0.35
    ) +
    geom_text(
        aes(label = star),
        size = 3.4,
        hjust = 0.5,
        vjust = 0.78,
        colour = outline,
        na.rm = TRUE
    ) +
    scale_fill_gradientn(
        name = config$labels$fill,
        colours = fill_cols,
        values = fill_values,
        limits = fill_lim,
        oob = scales::squish
    ) +
    scale_size_continuous(
        name = config$labels$size,
        range = c(1.6, 7.5)
    ) +
    guides(
        size = guide_legend(
            override.aes = list(
                fill = fill_cols[[5]],
                shape = 21,
                colour = outline,
                stroke = 0.35
            ),
            order = 1
        ),
        fill = guide_colorbar(order = 2)
    ) +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y
    ) +
    coord_cartesian(clip = "off") +
    theme_prism(base_size = 16, base_line_size = 0.85) +
    theme(
        panel.grid = element_blank(),
        plot.background = element_blank(),
        panel.background = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = outline, fill = NA, linewidth = 0.75),
        axis.line = element_blank(),
        axis.ticks = element_line(colour = outline, linewidth = 1.05),
        # Original template zeroed ticks (unit(0,"pt")); restore Prism-length ticks.
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(
            angle = 90, hjust = 1, vjust = 0.5,
            size = config$labels$x_size,
            margin = margin(t = 2.5)
        ),
        axis.text.y = element_text(
            size = config$labels$y_size,
            margin = margin(r = 2.5)
        ),
        strip.placement = "outside",
        strip.text.x = element_text(angle = 0, size = 15),
        strip.text.y = element_text(angle = 0, size = 15),
        strip.text.y.left = element_text(
            angle = 0, hjust = 1, vjust = 0.5,
            size = 15,
            margin = margin(r = 4)
        ),
        strip.text.y.right = element_text(angle = 0, size = 15),
        legend.key = element_blank(),
        panel.spacing.x = unit(0.7, "lines"),
        panel.spacing.y = unit(0.45, "lines")
    )

if (!is.null(row_f) || !is.null(col_f)) {
    facet_fml <- if (is.null(row_f)) {
        stats::as.formula(paste("~", col_f))
    } else if (is.null(col_f)) {
        stats::as.formula(paste(row_f, "~ ."))
    } else {
        stats::as.formula(paste(row_f, "~", col_f))
    }
    p <- p + facet_grid(
        facet_fml,
        scales = "free",
        space = "free",
        switch = "y"
    )
}

vline <- as.numeric(config$vline_intercepts)
vline <- vline[is.finite(vline)]
if (length(vline)) {
    p <- p + geom_vline(
        xintercept = vline,
        colour = outline,
        linewidth = 0.35,
        linetype = "22"
    )
}

if (!is.null(row_f)) {
    gt <- ggplotGrob(p)
    bar_col <- grid::segmentsGrob(
        x0 = unit(0.5, "npc"),
        x1 = unit(0.5, "npc"),
        y0 = unit(1.6, "mm"),
        y1 = unit(1, "npc") - unit(1.6, "mm"),
        gp = grid::gpar(col = outline, lwd = 2.6, lineend = "butt")
    )
    hits <- which(grepl("^strip-l", gt$layout$name))
    for (i in hits) {
        strip <- gt$grobs[[i]]
        if (!inherits(strip, "gtable")) {
            next
        }
        strip <- gtable::gtable_add_cols(strip, unit(1.6, "mm"), pos = -1)
        nc <- ncol(strip)
        strip <- gtable::gtable_add_grob(
            strip,
            bar_col,
            t = 1,
            b = nrow(strip),
            l = nc,
            r = nc,
            clip = "off",
            name = "group-bar"
        )
        strip <- gtable::gtable_add_cols(strip, unit(1.5, "mm"), pos = -1)
        gt$grobs[[i]] <- strip
    }
    p <- gt
}

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output, width = 12, height = 10)
