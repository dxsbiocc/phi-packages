#!/usr/bin/env Rscript

# Template-ID: bar-stack-donut
#
# Purpose:
#   Draw a horizontal absolute stacked bar chart of category frequencies
#   by a shared fill (e.g. gene × inheritance mode), with an inset
#   doughnut of overall fill composition (preview.png).
#
# Inputs:
#   Bars (example.tsv), one row per category-by-fill:
#     - gene, mode, count
#   Optional sidecar next to the input (donut.tsv):
#     - mode, percent   (or count; normalised for the doughnut)
#   If the sidecar is absent, doughnut values are the column sums of
#   the bar table (display aggregate of supplied counts).
#
# Output:
#   A PDF, PNG, or SVG stacked-bar + inset-doughnut figure.
#
# Dependencies:
#   ggplot2, readr, ggprism, scales, patchwork, grid
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   Edit CONFIG for columns, mode order/colours, donut sidecar name,
#   inset placement, labels, and size. Extra genes or modes are extra
#   rows/levels, not a new id. A plain stacked bar stays bar-stack; a
#   doughnut alone stays pie-doughnut; circos doughnut+bars stay
#   pie-doughnut-bar.
#
# Scientific assumptions:
#   Counts / percents are supplied display values. This script does not
#   call genotypes, assign inheritance, or select a gene list. Doughnut
#   slices without a sidecar are the sum of supplied bar counts per
#   mode (not an independent test).

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
        y = "gene",
        fill = "mode",
        value = "count"
    ),
    donut_file = "donut.tsv",
    donut_columns = list(
        fill = "mode",
        value = "percent"
    ),
    mode_order = c(
        "Autosomal Recessive",
        "Autosomal Recessive (de novo)",
        "Autosomal Dominant",
        "Autosomal Dominant (de novo)",
        "X-linked",
        "Non-Mendelian"
    ),
    palettes = list(fill = "Qualitative.Bold"),
    # 1-based picks into the fill palette (paper-like AR/AD/XL/NM).
    palette_index = c(5, 10, 2, 6, 3, 7),
    bar_height = 0.78,
    label_min = 1,
    italic_y = TRUE,
    donut_hole = 0.55,
    donut_label_min = 2.5,
    # patchwork inset_element npc box inside the bar panel.
    inset = list(left = 0.52, bottom = 0.12, right = 0.98, top = 0.72),
    labels = list(
        title = "63 Genes related to SNHL (Total diagnosed families, n=219)",
        x = "Frequency",
        y = "Gene name",
        donut_title = "Percentage of inheritance mode (n=220)",
        fill = NULL
    ),
    size = list(width = 12, height = 7.8)
)

load_packages(c("ggplot2", "readr", "ggprism", "scales", "patchwork", "grid"))

sidecar_path <- function(name) {
    if (is.null(name) || !nzchar(as.character(name)[[1]])) {
        return(NULL)
    }
    name <- as.character(name)[[1]]
    if (file.exists(name)) {
        return(normalizePath(name, mustWork = TRUE))
    }
    cand <- file.path(dirname(io$input), name)
    if (file.exists(cand)) {
        return(normalizePath(cand, mustWork = TRUE))
    }
    NULL
}

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- as.data.frame(read_table_auto(io$input), stringsAsFactors = FALSE)
require_columns(df, config$columns)
y_col <- config$columns$y
f_col <- config$columns$fill
v_col <- config$columns$value

df[[v_col]] <- as.numeric(df[[v_col]])
df <- df[is.finite(df[[v_col]]) & df[[v_col]] > 0, , drop = FALSE]
if (!nrow(df)) stop("No finite positive bar rows.", call. = FALSE)

mode_lv <- config$mode_order
mode_lv <- c(
    intersect(mode_lv, unique(as.character(df[[f_col]]))),
    setdiff(unique(as.character(df[[f_col]])), mode_lv)
)
df[[f_col]] <- factor(as.character(df[[f_col]]), levels = mode_lv)

tot <- tapply(df[[v_col]], as.character(df[[y_col]]), sum)
y_lv <- names(sort(tot, decreasing = TRUE))
df[[y_col]] <- factor(as.character(df[[y_col]]), levels = rev(y_lv))

tot_df <- data.frame(
    y = factor(names(tot), levels = levels(df[[y_col]])),
    total = as.numeric(tot),
    stringsAsFactors = FALSE
)
tot_df <- tot_df[tot_df$total >= config$label_min, , drop = FALSE]

n_m <- length(mode_lv)
pal <- palette_colors(config$palettes$fill)
idx <- as.integer(config$palette_index)
idx <- idx[is.finite(idx) & idx >= 1L & idx <= length(pal)]
fill_cols <- pal[idx]
if (length(fill_cols) < n_m) {
    fill_cols <- expand_palette(fill_cols, n_m)
} else {
    fill_cols <- fill_cols[seq_len(n_m)]
}
names(fill_cols) <- mode_lv

# Doughnut table: sidecar preferred; else sum of bar counts.
donut_path <- sidecar_path(config$donut_file)
if (!is.null(donut_path)) {
    dn <- as.data.frame(read_table_auto(donut_path), stringsAsFactors = FALSE)
    require_columns(dn, config$donut_columns)
    df_col <- config$donut_columns$fill
    dv_col <- config$donut_columns$value
    dn[[dv_col]] <- as.numeric(dn[[dv_col]])
    dn <- dn[is.finite(dn[[dv_col]]) & dn[[dv_col]] > 0, , drop = FALSE]
    dn$.mode <- as.character(dn[[df_col]])
    dn$.val <- dn[[dv_col]]
} else {
    agg <- tapply(df[[v_col]], as.character(df[[f_col]]), sum)
    dn <- data.frame(
        .mode = names(agg),
        .val = as.numeric(agg),
        stringsAsFactors = FALSE
    )
}
dn$.mode <- factor(dn$.mode, levels = mode_lv)
dn <- dn[!is.na(dn$.mode), , drop = FALSE]
dn <- dn[order(dn$.mode), , drop = FALSE]
dn$.pct <- 100 * dn$.val / sum(dn$.val)
dn$.ymax <- cumsum(dn$.pct)
dn$.ymin <- dn$.ymax - dn$.pct
dn$.mid <- (dn$.ymin + dn$.ymax) / 2
dn$.lab <- paste0(formatC(dn$.pct, format = "f", digits = 1), "%")
dn_lab <- dn[dn$.pct >= config$donut_label_min, , drop = FALSE]

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
fill_guide <- guide_legend(
    nrow = 1,
    byrow = TRUE,
    title.position = "top",
    override.aes = list(colour = NA)
)

p_bar <- ggplot(df, aes(
    x = .data[[v_col]],
    y = .data[[y_col]],
    fill = .data[[f_col]]
)) +
    geom_col(
        width = config$bar_height,
        colour = NA,
        position = position_stack(reverse = FALSE)
    ) +
    geom_label(
        data = tot_df,
        aes(x = total, y = y, label = total),
        inherit.aes = FALSE,
        hjust = -0.15,
        size = 2.6,
        linewidth = 0,
        fill = "white",
        label.padding = unit(0.1, "lines"),
        colour = "grey20"
    ) +
    scale_fill_manual(
        values = fill_cols,
        breaks = mode_lv,
        drop = FALSE,
        name = config$labels$fill,
        guide = fill_guide
    ) +
    scale_x_continuous(
        expand = expansion(mult = c(0, 0.12)),
        breaks = pretty_breaks(n = 6)
    ) +
    scale_y_discrete(expand = expansion(add = 0.35)) +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y
    ) +
    theme_prism(base_size = 11) +
    theme(
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.35),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(
            size = rel(0.72),
            face = if (isTRUE(config$italic_y)) "italic" else "plain",
            hjust = 1
        ),
        axis.text.x = element_text(size = rel(0.85)),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = rel(0.68)),
        legend.key.width = unit(0.85, "lines"),
        legend.key.height = unit(0.65, "lines"),
        plot.title = element_text(face = "bold", size = rel(1.0), hjust = 0),
        plot.background = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA),
        plot.margin = margin(8, 8, 4, 8)
    )

r_out <- 1
r_in <- config$donut_hole
p_donut <- ggplot(dn) +
    geom_rect(
        aes(xmin = r_in, xmax = r_out, ymin = .ymin, ymax = .ymax, fill = .mode),
        colour = "white",
        linewidth = 0.4
    ) +
    geom_text(
        data = dn_lab,
        aes(x = (r_in + r_out) / 2, y = .mid, label = .lab),
        size = 2.8,
        colour = "grey10"
    ) +
    coord_polar(theta = "y", start = 0, direction = -1) +
    xlim(0, r_out) +
    scale_fill_manual(values = fill_cols, guide = "none") +
    labs(title = config$labels$donut_title) +
    theme_void(base_size = 10) +
    theme(
        plot.title = element_text(
            face = "bold", size = rel(0.85), hjust = 0.5,
            margin = margin(b = 4)
        ),
        plot.background = element_rect(fill = "white", colour = NA),
        plot.margin = margin(2, 2, 2, 2)
    )

ins <- config$inset
p <- p_bar + inset_element(
    p_donut,
    left = ins$left,
    bottom = ins$bottom,
    right = ins$right,
    top = ins$top,
    align_to = "panel",
    on_top = TRUE,
    clip = FALSE
)

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(
    p,
    io$output,
    width = config$size$width,
    height = config$size$height
)
