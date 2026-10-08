#!/usr/bin/env Rscript

# Template-ID: line-ridge
#
# Purpose:
#   Draw stacked 1D density ridges (joyplot / Seurat RidgePlot): one
#   numeric axis, groups as overlapping ridgelines, optional facets.
#
# Inputs:
#   One row per observation. Default example:
#     - expression: numeric value (x)
#     - cluster: ridgeline group (y)
#     - marker: column facet
#     - lineage: optional row-group strip (Vγ9Vδ2 / Vδ1 / All cells)
#
# Output:
#   A PDF, PNG, or SVG ridgeline plot.
#
# Dependencies:
#   ggplot2, readr, ggprism, ggridges, gtable
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (columns, quantile, x limits,
#   labels). Facet by feature and left group bars are CONFIG, not a
#   second id. A single-panel ridge (no facet / no row_group) stays
#   on this id. Edit DATA PREPARATION to change group order. Edit
#   PLOT only when the glyph must change.
#
# Scientific assumptions:
#   Each row is one supplied observation of a numeric value in a
#   group. Ridges are a KDE display of those values. The vertical
#   tick is a display quantile of the same values. This script does
#   not cluster, call cell types, or pool observations into an
#   "All cells" group; that row is a supplied label.

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
        x = "expression",
        y = "cluster",
        fill = "cluster",
        facet = "marker",
        row_group = "lineage"
    ),
    quantile = 0.5,
    ridge_scale = 0.95,
    x_limits = c(0, 6),
    labels = list(
        title = "",
        x = "Expression",
        y = NULL
    )
)

load_packages(c("ggplot2", "readr", "ggprism", "ggridges", "gtable"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
need <- config$columns[c("x", "y")]
require_columns(df, need)

x_col <- config$columns$x
y_col <- config$columns$y
fill_col <- config$columns$fill
facet_col <- config$columns$facet
row_g <- config$columns$row_group

df[[x_col]] <- as.numeric(df[[x_col]])
ok <- is.finite(df[[x_col]])
if (any(!ok)) {
    message("Dropped ", sum(!ok), " row(s) with non-finite x.")
    df <- df[ok, , drop = FALSE]
}
if (!nrow(df)) {
    stop("No rows left to plot.", call. = FALSE)
}

if (is.null(fill_col) || !fill_col %in% names(df)) {
    fill_col <- y_col
} else if (!identical(fill_col, y_col)) {
    df[[fill_col]] <- factor(
        as.character(df[[fill_col]]),
        levels = unique(as.character(df[[fill_col]]))
    )
}

# File order of y is top-to-bottom (ggplot's first level sits at the bottom).
# Set last so a shared fill column cannot overwrite that order.
y_levels <- unique(as.character(df[[y_col]]))
df[[y_col]] <- factor(df[[y_col]], levels = rev(y_levels))

if (!is.null(facet_col) && facet_col %in% names(df)) {
    df[[facet_col]] <- factor(
        df[[facet_col]],
        levels = unique(as.character(df[[facet_col]]))
    )
} else {
    facet_col <- NULL
}

if (!is.null(row_g) && row_g %in% names(df)) {
    df[[row_g]] <- factor(
        df[[row_g]],
        levels = unique(as.character(df[[row_g]]))
    )
} else {
    row_g <- NULL
}

if (identical(fill_col, y_col)) {
    fill_lv <- y_levels
} else {
    fill_lv <- levels(df[[fill_col]])
}
alg <- palette_colors("Brand.Algolia")
safe <- palette_colors("Qualitative.Safe")
emma <- palette_colors("Brand.Emma")
# Paper-like cluster hues: blue, coral, green, purple, charcoal,
# peach, yellow, cyan. Extra groups continue Algolia.
ridge_seq <- c(
    alg[[8]], alg[[4]], alg[[2]], alg[[7]],
    emma[[1]], alg[[3]], safe[[3]], alg[[1]]
)
n_fill <- length(fill_lv)
if (n_fill <= length(ridge_seq)) {
    fill_vals <- ridge_seq[seq_len(n_fill)]
} else {
    fill_vals <- palette_colors("Brand.Algolia", n = n_fill)
    fill_vals[seq_along(ridge_seq)] <- ridge_seq
}
names(fill_vals) <- fill_lv
grey_hit <- grepl("all\\s*cells?", fill_lv, ignore.case = TRUE)
fill_vals[grey_hit] <- safe[[12]]

outline <- emma[[1]]
ink <- emma[[1]]

x_lim <- as.numeric(config$x_limits)
if (length(x_lim) != 2L || any(!is.finite(x_lim))) {
    x_lim <- range(df[[x_col]], na.rm = TRUE)
}
q_cut <- as.numeric(config$quantile)
q_cut <- q_cut[is.finite(q_cut) & q_cut > 0 & q_cut < 1]
use_q <- length(q_cut) > 0L
ridge_scale <- as.numeric(config$ridge_scale)
if (!is.finite(ridge_scale) || ridge_scale <= 0) {
    ridge_scale <- 0.95
}

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p <- ggplot(df, aes(
    x = .data[[x_col]],
    y = .data[[y_col]],
    fill = .data[[fill_col]]
)) +
    geom_density_ridges(
        scale = ridge_scale,
        rel_min_height = 0.01,
        from = x_lim[[1]],
        to = x_lim[[2]],
        quantile_lines = use_q,
        quantiles = if (use_q) q_cut else 0.5,
        colour = ink,
        vline_colour = ink,
        vline_width = 0.4,
        alpha = 0.92,
        linewidth = 0.25,
        show.legend = FALSE
    ) +
    scale_fill_manual(values = fill_vals, guide = "none") +
    scale_x_continuous(
        breaks = pretty(x_lim, n = 4),
        expand = expansion(mult = c(0.02, 0.04))
    ) +
    coord_cartesian(xlim = x_lim, clip = "off") +
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
        strip.background = element_blank(),
        panel.border = element_rect(colour = outline, fill = NA, linewidth = 0.55),
        axis.line = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.length.x = unit(2.5, "pt"),
        axis.text.y = element_text(hjust = 1),
        strip.placement = "outside",
        strip.text.x = element_text(angle = 0, face = "bold"),
        strip.text.y = element_blank(),
        strip.text.y.left = element_blank(),
        legend.key = element_blank(),
        panel.spacing.x = unit(0.45, "lines"),
        panel.spacing.y = unit(0.28, "lines")
    )

if (!is.null(row_g) || !is.null(facet_col)) {
    facet_fml <- if (is.null(row_g)) {
        stats::as.formula(paste("~", facet_col))
    } else if (is.null(facet_col)) {
        stats::as.formula(paste(row_g, "~ ."))
    } else {
        stats::as.formula(paste(row_g, "~", facet_col))
    }
    p <- p + facet_grid(
        facet_fml,
        scales = "free_y",
        space = "free_y"
    )
}

if (!is.null(row_g)) {
    gt <- ggplotGrob(p)
    bar_col <- grid::segmentsGrob(
        x0 = unit(0.5, "npc"),
        x1 = unit(0.5, "npc"),
        y0 = unit(1.6, "mm"),
        y1 = unit(1, "npc") - unit(1.6, "mm"),
        gp = grid::gpar(col = outline, lwd = 2.4, lineend = "butt")
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
save_ggplot(p, io$output, width = 9.8, height = 5.8)
