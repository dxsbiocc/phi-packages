#!/usr/bin/env Rscript

# Template-ID: line-ridge
#
# Purpose:
#   Draw stacked 1D density or histogram ridges (joyplot / Seurat
#   RidgePlot / FACS ratio ridges): one numeric axis, groups as
#   ridgelines, optional facets, optional per-group metric labels.
#
# Inputs:
#   1. Observations, one row per value. Default Seurat example:
#        - expression, cluster, marker, lineage
#      FACS example (example_facs.tsv):
#        - ratio, gene
#   2. Optional annotations table (gene, psi, mfi, ...), passed after
#      the observations and before the output file.
#
# Output:
#   A PDF, PNG, or SVG ridgeline plot.
#
# Dependencies:
#   ggplot2, readr, ggprism, ggridges, gtable
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#   RIDGE_EXAMPLE=facs Rscript plot.R example_facs.tsv annotations.tsv output.pdf
#   RIDGE_EXAMPLE=tcell Rscript plot.R example_tcell.tsv output.pdf
#
# Agent adaptation:
#   Edit only CONFIG (columns, quantile, x limits / log10, ridge
#   stat, fill palette, R label, facet_scales, y_labels). Facet by
#   feature and left group bars are CONFIG, not a second id.
#   Right-side PSI / MFI text comes from the optional annotations
#   table; R is a supplied display string. Edit DATA PREPARATION
#   to change group order.
#
# Scientific assumptions:
#   Each row is one supplied observation of a numeric value in a
#   group. Ridges are a KDE or histogram display of those values.
#   The vertical tick is a display quantile of the same values
#   when quantile is set. PSI, MFI, and R are supplied labels, not
#   computed here. This script does not cluster, call cell types,
#   or pool an "All cells" group; that row is a supplied label.

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

args <- commandArgs(trailingOnly = TRUE)
script <- "plot.R"
if (length(args) == 0L || args[[1]] %in% c("-h", "--help")) {
    cat(
        "Usage: Rscript plot.R <points> <output>\n",
        "       Rscript plot.R <points> <annotations> <output>\n",
        sep = ""
    )
    quit(save = "no", status = 0)
}
if (!length(args) %in% c(2L, 3L)) {
    stop(
        "Usage: Rscript plot.R <points> <output>\n",
        "       Rscript plot.R <points> <annotations> <output>",
        call. = FALSE
    )
}
if (!interactive() && unname(grDevices::dev.cur()) == 1L) {
    grDevices::pdf(nullfile())
}
points_path <- args[[1]]
ann_path <- if (length(args) == 3L) args[[2]] else NULL
out_path <- args[[length(args)]]

# -----------------------------------------------------------------------------
# CONFIG  (edit this block for a new dataset)
# -----------------------------------------------------------------------------
# Seurat multi-marker RidgePlot is the default. Presets:
#   RIDGE_EXAMPLE=facs   — GFP:mCherry FACS ridges + annotations.tsv
#   RIDGE_EXAMPLE=tcell  — T_MEM / T_RAC / T_EFF × gene facets
example_mode <- Sys.getenv("RIDGE_EXAMPLE", unset = "seurat")

config <- if (identical(example_mode, "facs")) {
    list(
        columns = list(
            x = "ratio",
            y = "gene",
            fill = "gene",
            facet = NULL,
            row_group = NULL,
            ann_y = "gene",
            ann_left = NULL,
            ann_right = c("psi", "mfi")
        ),
        quantile = NULL,
        ridge_scale = 0.92,
        ridge_stat = "density",
        bins = 45L,
        # Density is estimated on log10(x) when x_log10 is TRUE.
        bandwidth = 0.11,
        x_limits = c(1, 10000),
        x_log10 = TRUE,
        x_breaks = c(1, 10, 100, 1000, 10000),
        fill_palette = NULL,
        # Diverging.PiYG, magenta → green (skip near-white centre).
        fill_colours = c(
            CDC25A = "#c51b7d",
            FBXL14 = "#de77ae",
            USP18 = "#f1b6da",
            CKS1B = "#fde0ef",
            PARS2 = "#e6f5d0",
            TFEB = "#7fbc41",
            DNPEP = "#4d9221"
        ),
        r_label = "R = 0.98",
        r_label_x = 0.98,
        r_label_y = 0.98,
        annotation_size = 3.2,
        labels = list(
            title = "",
            x = "GFP:mCherry",
            y = NULL
        ),
        # Only the bottom x-axis line; no panel box.
        panel_border = FALSE,
        # Place group names on the first x tick (aligned with 10^0).
        y_at_first_break = TRUE,
        y_labels = NULL,
        facet_scales = "fixed",
        size = list(width = 7.4, height = 5.8)
    )
} else if (identical(example_mode, "tcell")) {
    list(
        columns = list(
            x = "expression",
            y = "cluster",
            fill = "cluster",
            facet = "marker",
            row_group = NULL,
            ann_y = NULL,
            ann_left = NULL,
            ann_right = NULL
        ),
        quantile = NULL,
        # Overlap so peaks sit above the next baseline (reference style).
        ridge_scale = 1.45,
        ridge_stat = "density",
        bins = 40L,
        bandwidth = 0.20,
        x_limits = NULL,
        x_log10 = FALSE,
        x_breaks = seq(-1, 6, by = 1),
        fill_palette = NULL,
        fill_colours = c(
            T_MEM = "#b84592",
            T_RAC = "#e6c547",
            T_EFF = "#4f83c4"
        ),
        r_label = NULL,
        r_label_x = 0.98,
        r_label_y = 0.98,
        annotation_size = 3.1,
        labels = list(
            title = "",
            x = "Expression Level",
            y = NULL
        ),
        panel_border = TRUE,
        y_at_first_break = FALSE,
        # Plotmath subscripts for T_MEM / T_RAC / T_EFF.
        y_labels = list(
            T_MEM = quote(T[MEM]),
            T_RAC = quote(T[RAC]),
            T_EFF = quote(T[EFF])
        ),
        facet_scales = "free_x",
        size = list(width = 8.6, height = 3.6)
    )
} else {
    list(
        columns = list(
            x = "expression",
            y = "cluster",
            fill = "cluster",
            facet = "marker",
            row_group = "lineage",
            ann_y = "gene",
            ann_left = NULL,
            ann_right = c("psi", "mfi")
        ),
        quantile = 0.5,
        ridge_scale = 0.95,
        ridge_stat = "density",
        bins = 40L,
        bandwidth = NULL,
        x_limits = c(0, 6),
        x_log10 = FALSE,
        x_breaks = NULL,
        fill_palette = NULL,
        fill_colours = NULL,
        r_label = NULL,
        r_label_x = 0.98,
        r_label_y = 0.98,
        annotation_size = 3.1,
        labels = list(
            title = "",
            x = "Expression",
            y = NULL
        ),
        panel_border = TRUE,
        y_at_first_break = FALSE,
        y_labels = NULL,
        facet_scales = "fixed",
        size = list(width = 9.8, height = 5.8)
    )
}

load_packages(c("ggplot2", "readr", "ggprism", "ggridges", "gtable", "scales"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- as.data.frame(read_table_auto(points_path), stringsAsFactors = FALSE)
need <- config$columns[c("x", "y")]
require_columns(df, need)

x_col <- config$columns$x
y_col <- config$columns$y
fill_col <- config$columns$fill
facet_col <- config$columns$facet
row_g <- config$columns$row_group

df[[x_col]] <- as.numeric(df[[x_col]])
ok <- is.finite(df[[x_col]])
if (isTRUE(config$x_log10)) {
    ok <- ok & df[[x_col]] > 0
}
if (any(!ok)) {
    message("Dropped ", sum(!ok), " row(s) with non-finite", if (isTRUE(config$x_log10)) " or non-positive", " x.")
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

ann <- NULL
if (!is.null(ann_path)) {
    ann <- as.data.frame(read_table_auto(ann_path), stringsAsFactors = FALSE)
    ann_y <- config$columns$ann_y
    if (is.null(ann_y) || !ann_y %in% names(ann)) {
        stop("Annotations table needs the ann_y column (default gene).", call. = FALSE)
    }
    ann[[ann_y]] <- as.character(ann[[ann_y]])
    right_cols <- config$columns$ann_right
    right_cols <- right_cols[!is.null(right_cols) & right_cols %in% names(ann)]
    if (!length(right_cols)) {
        stop("Annotations table needs at least one ann_right column.", call. = FALSE)
    }
    missing_y <- setdiff(y_levels, ann[[ann_y]])
    if (length(missing_y)) {
        stop(
            paste("Annotations missing groups:", paste(missing_y, collapse = ", ")),
            call. = FALSE
        )
    }
    ann <- ann[match(y_levels, ann[[ann_y]]), , drop = FALSE]
    ann$`.y` <- factor(ann[[ann_y]], levels = rev(y_levels))
    ann$`.label` <- apply(ann[right_cols], 1L, function(row) {
        parts <- vapply(seq_along(right_cols), function(i) {
            key <- toupper(right_cols[[i]])
            val <- row[[i]]
            if (is.numeric(val) || grepl("^[0-9.]+$", val)) {
                num <- as.numeric(val)
                if (is.finite(num) && abs(num - round(num)) < 1e-8) {
                    val <- as.character(as.integer(round(num)))
                } else if (is.finite(num)) {
                    val <- format(num, trim = TRUE, scientific = FALSE)
                }
            }
            paste0(key, " = ", val)
        }, character(1))
        paste(parts, collapse = "\n")
    })
}

n_fill <- length(fill_lv)
if (!is.null(config$fill_colours) && length(config$fill_colours)) {
    fill_vals <- config$fill_colours
    if (is.null(names(fill_vals))) {
        if (length(fill_vals) < n_fill) {
            stop("fill_colours is shorter than the number of groups.", call. = FALSE)
        }
        fill_vals <- setNames(fill_vals[seq_len(n_fill)], fill_lv)
    } else {
        miss <- setdiff(fill_lv, names(fill_vals))
        if (length(miss)) {
            stop(
                paste("fill_colours missing groups:", paste(miss, collapse = ", ")),
                call. = FALSE
            )
        }
        fill_vals <- fill_vals[fill_lv]
    }
} else if (!is.null(config$fill_palette) && nzchar(as.character(config$fill_palette))) {
    pal <- palette_colors(config$fill_palette, n = max(n_fill, 3L))
    fill_vals <- setNames(pal[seq_len(n_fill)], fill_lv)
} else {
    alg <- palette_colors("Brand.Algolia")
    safe <- palette_colors("Qualitative.Safe")
    emma <- palette_colors("Brand.Emma")
    ridge_seq <- c(
        alg[[8]], alg[[4]], alg[[2]], alg[[7]],
        emma[[1]], alg[[3]], safe[[3]], alg[[1]]
    )
    if (n_fill <= length(ridge_seq)) {
        fill_vals <- ridge_seq[seq_len(n_fill)]
    } else {
        fill_vals <- palette_colors("Brand.Algolia", n = n_fill)
        fill_vals[seq_along(ridge_seq)] <- ridge_seq
    }
    names(fill_vals) <- fill_lv
    grey_hit <- grepl("all\\s*cells?", fill_lv, ignore.case = TRUE)
    fill_vals[grey_hit] <- safe[[12]]
}

outline <- palette_colors("Brand.Emma")[[1]]
ink <- outline

x_lim <- config$x_limits
if (!is.null(x_lim)) {
    x_lim <- as.numeric(x_lim)
}
lock_xlim <- length(x_lim) == 2L && all(is.finite(x_lim))
if (!lock_xlim) {
    x_lim <- range(df[[x_col]], na.rm = TRUE)
    # Pad slightly so free_x panels still clear the density tails.
    pad <- 0.04 * diff(x_lim)
    if (!is.finite(pad) || pad <= 0) {
        pad <- 0.1
    }
    x_lim <- c(x_lim[[1]] - pad, x_lim[[2]] + pad)
}
q_cut <- as.numeric(config$quantile)
q_cut <- q_cut[is.finite(q_cut) & q_cut > 0 & q_cut < 1]
use_q <- length(q_cut) > 0L
ridge_scale <- as.numeric(config$ridge_scale)
if (!is.finite(ridge_scale) || ridge_scale <= 0) {
    ridge_scale <- 0.95
}
ridge_stat <- as.character(config$ridge_stat[[1]])
if (!ridge_stat %in% c("density", "binline")) {
    stop("ridge_stat must be \"density\" or \"binline\".", call. = FALSE)
}

# Log axes: estimate the density on log10(x). Using scale_x_log10 on
# raw ratios flattens linear-space KDEs into slabs.
plot_x <- x_col
plot_xlim <- x_lim
x_is_log <- isTRUE(config$x_log10)
if (x_is_log) {
    df$`.x_log` <- log10(df[[x_col]])
    plot_x <- ".x_log"
    plot_xlim <- log10(x_lim)
}

bw <- config$bandwidth
if (!is.null(bw)) {
    bw <- as.numeric(bw)[[1]]
    if (!is.finite(bw) || bw <= 0) {
        bw <- NULL
    }
}

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p <- ggplot(df, aes(
    x = .data[[plot_x]],
    y = .data[[y_col]],
    fill = .data[[fill_col]]
))

ridge_args <- list(
    scale = ridge_scale,
    rel_min_height = 0.001,
    colour = ink,
    alpha = 0.95,
    linewidth = 0.35,
    show.legend = FALSE
)
# When x limits are locked, clip the KDE to the axis window.
if (lock_xlim || x_is_log) {
    ridge_args$from <- plot_xlim[[1]]
    ridge_args$to <- plot_xlim[[2]]
}
if (!is.null(bw)) {
    ridge_args$bandwidth <- bw
}

if (identical(ridge_stat, "binline")) {
    p <- p +
        do.call(
            ggridges::geom_density_ridges,
            c(
                list(
                    stat = "binline",
                    bins = as.integer(config$bins)
                ),
                ridge_args[setdiff(names(ridge_args), "bandwidth")]
            )
        )
} else {
    p <- p +
        do.call(
            ggridges::geom_density_ridges,
            c(
                list(
                    quantile_lines = use_q,
                    quantiles = if (use_q) q_cut else 0.5,
                    vline_colour = ink,
                    vline_width = 0.4
                ),
                ridge_args
            )
        )
}

# Baseline across the full panel width for every ridge (-Inf/Inf so
# free_x facets still span left-to-right).
base_df <- data.frame(
    y = factor(y_levels, levels = rev(y_levels)),
    stringsAsFactors = FALSE
)
p <- p +
    geom_segment(
        data = base_df,
        aes(y = y, yend = y),
        x = -Inf,
        xend = Inf,
        inherit.aes = FALSE,
        colour = ink,
        linewidth = 0.4,
        lineend = "butt"
    )

p <- p +
    scale_fill_manual(values = fill_vals, guide = "none")

brks <- config$x_breaks
if (x_is_log) {
    if (is.null(brks) || !length(brks)) {
        brks <- 10^seq(
            floor(plot_xlim[[1]]),
            ceiling(plot_xlim[[2]]),
            by = 1
        )
    }
    brks <- as.numeric(brks)
    p <- p +
        scale_x_continuous(
            breaks = log10(brks),
            labels = scales::math_format(10^.x)(log10(brks)),
            expand = expansion(mult = c(0, 0.10))
        )
    first_break <- log10(brks[[1]])
} else {
    if (is.null(brks) || !length(brks)) {
        brks <- pretty(plot_xlim, n = 4)
    }
    brks <- as.numeric(brks)
    p <- p +
        scale_x_continuous(
            breaks = brks,
            expand = expansion(mult = c(0.02, 0.04))
        )
    first_break <- brks[[1]]
}

# Bottom flush with the x-axis; small top pad for the uppermost peak.
y_expand <- expansion(mult = c(0, 0), add = c(0, 0.55))
y_lab_map <- config$y_labels
if (!is.null(y_lab_map) && length(y_lab_map)) {
    # Named plotmath / character map; missing keys fall back to the level.
    y_lab_fun <- function(breaks) {
        lapply(breaks, function(b) {
            if (!is.null(y_lab_map[[b]])) y_lab_map[[b]] else b
        })
    }
} else {
    y_lab_fun <- waiver()
}

y_inside <- isTRUE(config$y_at_first_break)
if (y_inside) {
    p <- p +
        scale_y_discrete(
            breaks = levels(df[[y_col]]),
            labels = NULL,
            expand = y_expand
        )
    y_lab <- data.frame(
        lab_x = first_break,
        lab_y = factor(y_levels, levels = rev(y_levels)),
        lab = y_levels,
        stringsAsFactors = FALSE
    )
    p <- p +
        geom_text(
            data = y_lab,
            aes(x = .data$lab_x, y = .data$lab_y, label = .data$lab),
            inherit.aes = FALSE,
            hjust = 0,
            # Sit just above the full-width baseline.
            vjust = 0,
            nudge_y = 0.08,
            size = as.numeric(config$annotation_size) + 0.5,
            colour = ink,
            fontface = "bold"
        )
} else {
    p <- p +
        scale_y_discrete(labels = y_lab_fun, expand = y_expand)
}

base_theme <- if (y_inside) {
    theme_classic(base_size = 12)
} else {
    theme_prism()
}

if (lock_xlim || x_is_log) {
    p <- p + coord_cartesian(xlim = plot_xlim, clip = "off")
} else {
    p <- p + coord_cartesian(clip = "off")
}

p <- p +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y
    ) +
    base_theme +
    theme(
        panel.grid = element_blank(),
        plot.background = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA),
        legend.background = element_blank(),
        strip.background = element_blank(),
        panel.border = if (isTRUE(config$panel_border)) {
            element_rect(colour = outline, fill = NA, linewidth = 0.55)
        } else {
            element_blank()
        },
        axis.line.x = if (isTRUE(config$panel_border)) {
            element_blank()
        } else {
            element_line(colour = ink, linewidth = 0.55)
        },
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.length.y = unit(0, "pt"),
        axis.ticks.length.x = unit(2.5, "pt"),
        axis.text.y = if (y_inside) element_blank() else element_text(hjust = 1, colour = ink),
        axis.text.y.left = if (y_inside) element_blank() else element_text(hjust = 1, colour = ink),
        axis.text.x = element_text(colour = ink),
        axis.title = element_text(colour = ink),
        axis.title.y = element_blank(),
        strip.placement = "outside",
        strip.text.x = element_text(angle = 0, face = "bold", colour = ink),
        strip.text.y = element_blank(),
        strip.text.y.left = element_blank(),
        legend.key = element_blank(),
        panel.spacing.x = unit(0.55, "lines"),
        panel.spacing.y = unit(0.28, "lines"),
        plot.margin = if (y_inside) margin(8, 56, 8, 8) else margin(5.5, 5.5, 5.5, 5.5)
    )

if (!is.null(ann)) {
    # Place labels just past the right of the data window.
    x_lab <- plot_xlim[[2]] + 0.06 * diff(plot_xlim)
    ann$`.x` <- x_lab
    p <- p +
        geom_text(
            data = ann,
            aes(x = .data$`.x`, y = .data$`.y`, label = .data$`.label`),
            inherit.aes = FALSE,
            hjust = 0,
            vjust = 0.5,
            size = as.numeric(config$annotation_size),
            colour = ink,
            lineheight = 0.92
        ) +
        theme(plot.margin = margin(8, 56, 8, 8))
}

if (!is.null(config$r_label) && nzchar(as.character(config$r_label))) {
    p <- p +
        annotate(
            "text",
            x = plot_xlim[[1]] + as.numeric(config$r_label_x) * diff(plot_xlim),
            y = Inf,
            label = as.character(config$r_label),
            hjust = 1,
            vjust = 1.4,
            size = as.numeric(config$annotation_size) + 0.4,
            colour = ink
        )
}

if (!is.null(row_g) || !is.null(facet_col)) {
    facet_fml <- if (is.null(row_g)) {
        stats::as.formula(paste("~", facet_col))
    } else if (is.null(facet_col)) {
        stats::as.formula(paste(row_g, "~ ."))
    } else {
        stats::as.formula(paste(row_g, "~", facet_col))
    }
    facet_scales <- as.character(
        if (is.null(config$facet_scales)) "fixed" else config$facet_scales
    )[[1]]
    if (!facet_scales %in% c("fixed", "free", "free_x", "free_y")) {
        facet_scales <- "fixed"
    }
    # Row groups still need free_y / free_y space so lineage blocks size
    # to their cluster counts; free_x overlays for per-gene axes.
    if (!is.null(row_g)) {
        if (identical(facet_scales, "fixed")) {
            facet_scales <- "free_y"
        } else if (identical(facet_scales, "free_x")) {
            facet_scales <- "free"
        }
        p <- p + facet_grid(
            facet_fml,
            scales = facet_scales,
            space = "free_y"
        )
    } else {
        p <- p + facet_grid(
            facet_fml,
            scales = facet_scales,
            space = "fixed"
        )
    }
}

if (!is.null(row_g) || y_inside) {
    gt <- ggplotGrob(p)
    if (!is.null(row_g)) {
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
    }
    if (y_inside) {
        # Drop leftover left-axis / ylab columns so names sit on the
        # first x tick inside the panel, not in an outer margin.
        drop_cols <- unique(unlist(lapply(
            c("axis-l", "ylab-l"),
            function(nm) {
                hit <- gt$layout$l[gt$layout$name == nm]
                if (!length(hit)) integer(0) else hit
            }
        )))
        for (col in sort(unique(drop_cols), decreasing = TRUE)) {
            gt$widths[[col]] <- unit(0, "pt")
        }
    }
    p <- gt
}

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(
    p,
    out_path,
    width = config$size$width,
    height = config$size$height
)
