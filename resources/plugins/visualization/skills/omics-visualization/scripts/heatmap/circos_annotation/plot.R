#!/usr/bin/env Rscript

# Template-ID: heatmap-circos-annotation
#
# Purpose:
#   Circular cohort overview: one equal-angle slice per unit (sample /
#   tissue block), concentric rings of independent categorical
#   attributes (e.g. cohort, assay, tumor type, cancer type), an
#   optional outer radial bar per unit (e.g. number of serial
#   sections), and the totals in the centre. An optional fan-shaped
#   ring key explains every ring.
#
# Inputs:
#   One row per unit (tsv/csv). Default example:
#     - cancer, tumor, assay, cohort: categorical ring columns
#     - n_sections: optional numeric bar per unit
#
# Output:
#   A PDF, PNG, or SVG circular annotation figure.
#
# Dependencies:
#   ggplot2, readr, gground, patchwork
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   Edit only CONFIG: rings (column -> ring title, inner to outer),
#   label_ring, bar column, colours, ring_key, sizes.
#
# Scientific assumptions:
#   Rings are independent supplied attributes, not a hierarchy, so ring
#   order carries no meaning beyond layout. Units are sorted by
#   label_ring then the other rings, and adjacent equal values merge
#   into one arc. Labels such as "BRCA 54(31)" and the centre "122 (98)"
#   are tallies of the supplied rows (sum of bar, number of units);
#   nothing is estimated.

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
    # Inner to outer; names are columns, values are ring titles.
    rings = c(
        cohort = "Cohort",
        assay = "Assay",
        tumor = "Tumor type",
        cancer = "Cancer type"
    ),
    label_ring = "cancer",      # thicker ring with "<value>\n<bar sum>(<n>)" labels
    bar = "n_sections",         # NULL = no outer bars
    bar_title = "n of sections",
    unit_title = "n of tissues",
    # Named hex per ring value; unlisted values take the palette.
    colours = list(
        cancer = c(BRCA = "#d7689d", CHOL = "#f3d337", CRC = "#a0c443",
            PDAC = "#66b5ae", RCC = "#5da4dc", UCEC = "#40769e"),
        tumor = c(Metastasis = "#f26b21", Primary = "#68b8be"),
        assay = c(OCT = "#168de2", FFPE = "#ffa930"),
        cohort = c(Distinct = "#00457c", Diffuse = "#e8d4a2")
    ),
    bar_colours = c("#eff3ff", "#bdd7e7", "#6baed6", "#3182bd", "#08519c"),
    palettes = list(ring = "Qualitative.Safe"),
    ring_width = 0.9,
    label_ring_width = 1.8,
    bar_gap = 0.2,
    bar_scale = 1,
    corner = grid::unit(3, "pt"),
    label_size = 3.6,
    centre_size = 5,
    ring_key = TRUE,            # fan-shaped key beside the circle
    widths = c(1.45, 1)
)

load_packages(c("ggplot2", "readr", "gground", "patchwork"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- as.data.frame(read_table_auto(io$input), stringsAsFactors = FALSE)
ring_cols <- names(config$rings)
if (!config$label_ring %in% ring_cols) {
    stop("label_ring must be one of the ring columns.", call. = FALSE)
}
has_bar <- !is.null(config$bar) && config$bar %in% names(df)
require_columns(df, as.list(c(ring_cols, if (has_bar) config$bar)))
for (cl in ring_cols) {
    df[[cl]] <- as.character(df[[cl]])
}
df <- df[stats::complete.cases(df[ring_cols]), , drop = FALSE]
if (has_bar) {
    df$.bar <- as.numeric(df[[config$bar]])
    df$.bar[!is.finite(df$.bar)] <- 0
} else {
    df$.bar <- 0
}
if (!nrow(df)) {
    stop("No complete rows to plot.", call. = FALSE)
}

# Sort: label ring first, then the other rings outer -> inner, then bar.
sort_cols <- c(config$label_ring, rev(setdiff(ring_cols, config$label_ring)))
ord <- do.call(order, c(lapply(sort_cols, function(cl) df[[cl]]), list(-df$.bar)))
df <- df[ord, , drop = FALSE]
n_unit <- nrow(df)

# Colours per ring value (supplied first, palette for the rest).
ring_fill <- lapply(ring_cols, function(cl) {
    lv <- unique(df[[cl]])
    fixed <- config$colours[[cl]]
    fixed <- if (is.null(fixed)) character(0) else fixed[names(fixed) %in% lv]
    miss <- setdiff(lv, names(fixed))
    if (length(miss)) {
        extra <- palette_colors(config$palettes$ring, n = length(miss))
        names(extra) <- miss
        fixed <- c(fixed, extra)
    }
    fixed[lv]
})
names(ring_fill) <- ring_cols
bar_levels <- sort(unique(df$.bar[df$.bar > 0]))
bar_fill <- expand_palette(config$bar_colours, max(1, length(bar_levels)))
names(bar_fill) <- bar_levels

# Ring x positions (radius): inner rings, then the label ring, then bars.
x0 <- numeric(length(ring_cols))
x1 <- numeric(length(ring_cols))
pos <- 1
for (k in seq_along(ring_cols)) {
    w <- if (ring_cols[[k]] == config$label_ring) config$label_ring_width else config$ring_width
    x0[[k]] <- pos
    x1[[k]] <- pos + w
    pos <- pos + w + (1 - config$ring_width)
}
bar_x0 <- max(x1) + config$bar_gap + 0.1

runs_of <- function(v) {
    r <- rle(v)
    ymax <- cumsum(r$lengths)
    data.frame(ymin = ymax - r$lengths, ymax = ymax, value = r$values,
        stringsAsFactors = FALSE)
}
rects <- do.call(rbind, lapply(seq_along(ring_cols), function(k) {
    cl <- ring_cols[[k]]
    d <- runs_of(df[[cl]])
    d$xmin <- x0[[k]]
    d$xmax <- x1[[k]]
    d$fill <- unname(ring_fill[[cl]][d$value])
    d$ring <- cl
    d
}))
bars <- NULL
if (has_bar) {
    b <- runs_of(df$.bar)
    b <- b[b$value > 0, , drop = FALSE]
    b$xmin <- bar_x0
    b$xmax <- bar_x0 + b$value * config$bar_scale
    b$fill <- unname(bar_fill[as.character(b$value)])
    bars <- b
}

tangent <- function(y, total) {
    a <- -360 * y / total
    a <- ((a + 180) %% 360) - 180
    ifelse(a > 90 | a < -90, a + 180, a)
}
lab_ring <- rects[rects$ring == config$label_ring, , drop = FALSE]
# Merge runs by value so each category gets one label.
lab <- do.call(rbind, lapply(unique(lab_ring$value), function(v) {
    s <- lab_ring[lab_ring$value == v, , drop = FALSE]
    rows <- df[df[[config$label_ring]] == v, , drop = FALSE]
    n_txt <- if (has_bar) sprintf("%g(%d)", sum(rows$.bar), nrow(rows)) else sprintf("%d", nrow(rows))
    data.frame(
        x = (s$xmin[[1]] + s$xmax[[1]]) / 2,
        y = (min(s$ymin) + max(s$ymax)) / 2,
        label = paste0(v, "\n", n_txt),
        stringsAsFactors = FALSE
    )
}))
lab$angle <- tangent(lab$y, n_unit)
centre <- if (has_bar) sprintf("%g\n(%d)", sum(df$.bar), n_unit) else sprintf("%d", n_unit)
x_lim <- c(-1, if (has_bar) bar_x0 + max(df$.bar) * config$bar_scale + 0.3 else max(x1) + 0.5)

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
# One layer per fill colour: geom_round_rect scrambles fills when one
# layer mixes many colours under coord_polar.
draw_rects <- function(p, d) {
    for (f in unique(d$fill)) {
        p <- p + gground::geom_round_rect(
            data = d[d$fill == f, , drop = FALSE],
            ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = f, radius = config$corner, colour = NA
        )
    }
    p
}
p_main <- ggplot2::ggplot()
p_main <- draw_rects(p_main, rects)
if (!is.null(bars)) p_main <- draw_rects(p_main, bars)
p_main <- p_main +
    ggplot2::geom_text(
        data = lab, ggplot2::aes(x = x, y = y, label = label, angle = angle),
        colour = "white", fontface = "bold", size = config$label_size, lineheight = 0.85
    ) +
    ggplot2::annotate("text", x = -1, y = 0, label = centre,
        size = config$centre_size, lineheight = 0.85) +
    ggplot2::scale_fill_identity() +
    ggplot2::coord_polar(theta = "y", clip = "off") +
    ggplot2::scale_x_continuous(limits = x_lim) +
    ggplot2::scale_y_continuous(limits = c(0, n_unit), expand = c(0, 0)) +
    ggplot2::theme_void() +
    ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white", colour = NA))

p <- p_main
if (isTRUE(config$ring_key)) {
    # Fan key: each ring spans a half circle, values split equally.
    half <- 6
    key <- do.call(rbind, lapply(seq_along(ring_cols), function(k) {
        cl <- ring_cols[[k]]
        lv <- names(ring_fill[[cl]])
        w <- half / length(lv)
        data.frame(
            xmin = x0[[k]], xmax = x1[[k]],
            ymin = (seq_along(lv) - 1) * w + 0.04, ymax = seq_along(lv) * w - 0.04,
            fill = unname(ring_fill[[cl]]), value = lv, stringsAsFactors = FALSE
        )
    }))
    key$x <- (key$xmin + key$xmax) / 2
    key$y <- (key$ymin + key$ymax) / 2
    key$angle <- tangent(key$y, 2 * half)
    titles <- data.frame(
        x = c(-1, (x0 + x1) / 2),
        label = c(sprintf("%s\n(%s)", config$bar_title, config$unit_title), unname(config$rings)),
        stringsAsFactors = FALSE
    )
    p_key <- ggplot2::ggplot()
    p_key <- draw_rects(p_key, key)
    key_bar_scale <- 0.5 * config$bar_scale
    key_xmax <- max(x1) + 0.5
    if (has_bar && length(bar_levels)) {
        key_xmax <- bar_x0 + max(bar_levels) * key_bar_scale + 0.3
        kb <- data.frame(
            xmin = bar_x0,
            xmax = bar_x0 + rev(bar_levels) * key_bar_scale,
            ymin = (seq_along(bar_levels) - 1) * 0.55 + 0.05,
            ymax = (seq_along(bar_levels) - 1) * 0.55 + 0.45,
            fill = unname(bar_fill[as.character(rev(bar_levels))]),
            value = rev(bar_levels)
        )
        p_key <- draw_rects(p_key, kb) +
            ggplot2::geom_text(
                data = kb,
                ggplot2::aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = value),
                size = 3.2, colour = ifelse(kb$value >= stats::median(bar_levels), "white", "grey15")
            )
        titles <- rbind(titles, data.frame(x = bar_x0 + 0.3, label = config$bar_title))
    }
    p_key <- p_key +
        ggplot2::geom_text(
            data = key, ggplot2::aes(x = x, y = y, label = value, angle = angle),
            colour = "white", fontface = "bold", size = 3
        ) +
        ggplot2::geom_text(
            data = titles, ggplot2::aes(x = x, y = 0, label = label),
            hjust = 1.06, size = 3.6, lineheight = 0.85
        ) +
        ggplot2::annotate("text", x = -1, y = 0.9, label = centre,
            size = config$centre_size - 1, lineheight = 0.85, hjust = 0) +
        ggplot2::scale_fill_identity() +
        ggplot2::coord_polar(theta = "y", clip = "off") +
        ggplot2::scale_x_continuous(limits = c(-1, key_xmax)) +
        ggplot2::scale_y_continuous(limits = c(0, 2 * half), expand = c(0, 0)) +
        ggplot2::theme_void() +
        ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white", colour = NA))
    p <- patchwork::wrap_plots(p_main, p_key, nrow = 1, widths = config$widths)
}

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output, width = if (isTRUE(config$ring_key)) 12.5 else 7.5, height = 7.2)
