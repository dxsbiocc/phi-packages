#!/usr/bin/env Rscript

# Template-ID: bar-stack
#
# Purpose:
#   Draw stacked bars of a component value within each category, either
#   normalised to 100% (default) or as absolute magnitudes, with optional
#   translucent bands linking the same component between adjacent bars.
#
# Inputs:
#   A table with one row per category-by-component. Default example:
#     - day: category (x-axis, file order)
#     - value: numeric value
#     - category: component (stack layer, file order = bottom to top)
#
# Output:
#   A PDF, PNG, or SVG stacked bar chart.
#
# Dependencies:
#   ggplot2, readr, ggprism
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   Edit only CONFIG: columns, labels, scale ("percent" or "absolute"),
#   links (bands between bars), bar width, label rules, palette.
#   Edit PLOT only when the chart geometry or styling must change.
#
# Scientific assumptions:
#   Each row is one category-by-component pair with a comparable value.
#   With scale = "percent" each bar is divided by its own total, so bars
#   show composition, not size. Link bands only connect the same
#   component in neighbouring bars to guide the eye; they are not flows
#   between components (that would be a Sankey). Values are already
#   summarised, not raw replicates.

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
        x = "day",
        y = "value",
        fill = "category"
    ),
    labels = list(
        title = "",
        x = "",
        y = NULL                # NULL = "Proportion" or "Value" by scale
    ),
    scale = "percent",          # "percent" (each bar sums to 1) or "absolute"
    links = TRUE,               # translucent bands joining the same component
    link_alpha = 0.28,
    bar_width = 0.6,
    value_labels = TRUE,
    label_min = 0.03,           # hide labels on segments smaller than this share
    label_digits = 1,
    palettes = list(fill = "Qualitative.Prism"),
    # Palette positions in use order, so adjacent layers contrast.
    palette_order = c(2, 5, 9, 7, 3, 6, 1, 8, 4, 10, 11),
    legend_position = "bottom"
)

load_packages(c("ggplot2", "readr", "ggprism"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- as.data.frame(read_table_auto(io$input), stringsAsFactors = FALSE)
require_columns(df, config$columns)
x_col <- config$columns$x
y_col <- config$columns$y
fill_col <- config$columns$fill
df$.x <- factor(df[[x_col]], levels = unique(df[[x_col]]))
df$.f <- factor(df[[fill_col]], levels = unique(df[[fill_col]]))
df$.v <- as.numeric(df[[y_col]])
df <- df[is.finite(df$.v) & df$.v >= 0, , drop = FALSE]
if (!nrow(df)) {
    stop("No finite non-negative values to plot.", call. = FALSE)
}
if (anyDuplicated(paste(df$.x, df$.f, sep = "\r"))) {
    stop("Each category-by-component pair must appear once.", call. = FALSE)
}
scale_mode <- match.arg(config$scale, c("percent", "absolute"))

df <- df[order(df$.x, df$.f), , drop = FALSE]
df$.xi <- as.integer(df$.x)
df$share <- stats::ave(df$.v, df$.xi, FUN = function(v) v / sum(v))
df$h <- if (identical(scale_mode, "percent")) df$share else df$.v
df$ymax <- stats::ave(df$h, df$.xi, FUN = cumsum)
df$ymin <- df$ymax - df$h

w <- config$bar_width
f_levels <- levels(df$.f)
pal <- palette_colors(config$palettes$fill)
ord <- config$palette_order[config$palette_order <= length(pal)]
pal <- c(pal[ord], pal[setdiff(seq_along(pal), ord)])
fill_cols <- expand_palette(pal, length(f_levels))
names(fill_cols) <- f_levels

# Link bands: right edge of bar i to left edge of bar i + 1, per component.
poly <- NULL
if (isTRUE(config$links) && nlevels(df$.x) > 1L) {
    pieces <- list()
    for (lv in f_levels) {
        d <- df[df$.f == lv, , drop = FALSE]
        d <- d[order(d$.xi), , drop = FALSE]
        if (nrow(d) < 2L) next
        for (k in seq_len(nrow(d) - 1L)) {
            a <- d[k, ]
            b <- d[k + 1L, ]
            if (b$.xi != a$.xi + 1L) next
            pieces[[length(pieces) + 1L]] <- data.frame(
                px = c(a$.xi + w / 2, b$.xi - w / 2, b$.xi - w / 2, a$.xi + w / 2),
                py = c(a$ymin, b$ymin, b$ymax, a$ymax),
                id = paste(lv, a$.xi),
                fill = fill_cols[[lv]],
                stringsAsFactors = FALSE
            )
        }
    }
    if (length(pieces)) poly <- do.call(rbind, pieces)
}

text_on <- function(hex) {
    rgb <- grDevices::col2rgb(hex) / 255
    lum <- 0.2126 * rgb[1, ] + 0.7152 * rgb[2, ] + 0.0722 * rgb[3, ]
    ifelse(lum > 0.6, "grey15", "white")
}
lab <- NULL
if (isTRUE(config$value_labels)) {
    lab <- df[df$share >= config$label_min, , drop = FALSE]
    lab$label <- if (identical(scale_mode, "percent")) {
        paste0(formatC(100 * lab$share, format = "f", digits = config$label_digits), "%")
    } else {
        format(lab$.v, big.mark = ",", trim = TRUE)
    }
    lab$col <- text_on(fill_cols[as.character(lab$.f)])
}

y_lab <- config$labels$y
if (is.null(y_lab)) {
    y_lab <- if (identical(scale_mode, "percent")) "Proportion" else "Value"
}

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p <- ggplot2::ggplot()
if (!is.null(poly)) {
    p <- p + ggplot2::geom_polygon(
        data = poly, ggplot2::aes(x = px, y = py, group = id),
        fill = poly$fill, alpha = config$link_alpha, colour = NA
    )
}
p <- p +
    ggplot2::geom_rect(
        data = df,
        ggplot2::aes(xmin = .xi - w / 2, xmax = .xi + w / 2, ymin = ymin, ymax = ymax,
            fill = .f),
        colour = NA
    )
if (!is.null(lab)) {
    p <- p + ggplot2::geom_text(
        data = lab,
        ggplot2::aes(x = .xi, y = (ymin + ymax) / 2, label = label),
        colour = lab$col, size = 3.3
    )
}
p <- p +
    ggplot2::scale_fill_manual(values = fill_cols, name = NULL) +
    ggplot2::scale_x_continuous(
        breaks = seq_len(nlevels(df$.x)), labels = levels(df$.x),
        expand = ggplot2::expansion(add = 0.5)
    ) +
    ggplot2::scale_y_continuous(
        expand = ggplot2::expansion(mult = c(0, 0.02)),
        labels = if (identical(scale_mode, "percent")) {
            function(v) paste0(round(100 * v), "%")
        } else {
            ggplot2::waiver()
        }
    ) +
    ggplot2::labs(title = config$labels$title, x = config$labels$x, y = y_lab) +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1, byrow = TRUE)) +
    ggprism::theme_prism(base_size = 12) +
    ggplot2::theme(
        legend.position = config$legend_position,
        legend.text = ggplot2::element_text(size = 10, face = "plain"),
        panel.grid.major.y = ggplot2::element_line(colour = "grey90", linewidth = 0.4),
        plot.background = ggplot2::element_rect(fill = "white", colour = NA),
        plot.margin = ggplot2::margin(10, 12, 8, 8)
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output, width = max(6.5, 0.95 * nlevels(df$.x) + 1.8), height = 6)
