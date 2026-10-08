#!/usr/bin/env Rscript

# Template-ID: heatmap-mutation-energy
#
# Purpose:
#   Draw a saturation-mutation energy heatmap (site × amino acid) with
#   a dual-y overlay of mean ΔΔG per site.
#
# Inputs:
#   One row per residue. Default example:
#     - site: residue label
#     - remaining columns: one numeric ΔΔG per mutant amino acid
#
# Output:
#   A PDF, PNG, or SVG heatmap with a mean-energy line.
#
# Dependencies:
#   ggplot2, readr, ggprism, scales, aplot, patchwork
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (site column, AA order, labels,
#   fill limits). Edit DATA PREPARATION to reshape a long table.
#   Edit PLOT only when the overlay geometry must change.
#   This is one template, not a heatmap plus a separate line chart.
#
# Scientific assumptions:
#   Each cell is a supplied ΔΔG for one substitution. This script does
#   not dock, fold, or estimate energies. The overlay line is the mean
#   of the supplied cells in that site (a display aggregate, not a new
#   assay). Do not interpolate tiles; one cell is one mutation.

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
        site = "site"
    ),
    aa_order = c(
        "A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
        "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"
    ),
    fill_limits = c(-1.5, 3),
    labels = list(
        title = "",
        x = NULL,
        y = "Mutation",
        y_right = "Average mutation energy\n(binding, kcal/mol)",
        fill = "\u0394\u0394G (kcal/mol)"
    )
)

load_packages(c("ggplot2", "readr", "ggprism", "scales", "aplot", "patchwork"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
require_columns(df, config$columns)
site_col <- config$columns$site
aa_cols <- setdiff(names(df), site_col)
if (!length(aa_cols)) {
    stop("Need at least one amino-acid ΔΔG column besides site.", call. = FALSE)
}
wanted <- as.character(config$aa_order)
aa_cols <- c(intersect(wanted, aa_cols), setdiff(aa_cols, wanted))
for (col in aa_cols) {
    df[[col]] <- as.numeric(df[[col]])
}
ok <- !is.na(df[[site_col]]) & nzchar(as.character(df[[site_col]]))
if (any(!ok)) {
    message("Dropped ", sum(!ok), " row(s) with missing site labels.")
    df <- df[ok, , drop = FALSE]
}
if (!nrow(df)) {
    stop("No sites left to plot.", call. = FALSE)
}
df[[site_col]] <- factor(df[[site_col]], levels = unique(as.character(df[[site_col]])))
n_site <- nrow(df)
n_aa <- length(aa_cols)

heat <- data.frame(
    site = rep(df[[site_col]], times = n_aa),
    aa = factor(rep(aa_cols, each = n_site), levels = aa_cols),
    ddg = as.numeric(as.matrix(df[aa_cols])),
    stringsAsFactors = FALSE
)
heat$x <- as.numeric(heat$site)
heat$y <- as.numeric(heat$aa)

avg <- data.frame(
    site = df[[site_col]],
    x = seq_len(n_site),
    mean_ddg = rowMeans(as.matrix(df[aa_cols]), na.rm = TRUE),
    stringsAsFactors = FALSE
)
avg <- avg[is.finite(avg$mean_ddg), , drop = FALSE]
if (!nrow(avg)) {
    stop("No finite site means to overlay.", call. = FALSE)
}

fill_lim <- as.numeric(config$fill_limits)
if (length(fill_lim) != 2L || any(!is.finite(fill_lim))) {
    fill_lim <- range(heat$ddg, na.rm = TRUE)
}
ddg_low <- min(c(avg$mean_ddg, 0), na.rm = TRUE)
ddg_hig <- max(c(avg$mean_ddg, 0), na.rm = TRUE)
if (!is.finite(ddg_low) || !is.finite(ddg_hig) || ddg_low == ddg_hig) {
    ddg_low <- ddg_low - 0.5
    ddg_hig <- ddg_hig + 0.5
}
linmap <- function(x, from, to) {
    to[[1]] + (x - from[[1]]) / (from[[2]] - from[[1]]) * (to[[2]] - to[[1]])
}
avg$y_line <- linmap(avg$mean_ddg, c(ddg_low, ddg_hig), c(1, n_aa))
avg$hot <- avg$mean_ddg > 0
y_zero <- linmap(0, c(ddg_low, ddg_hig), c(1, n_aa))
to_ddg <- function(y) {
    linmap(y, c(1, n_aa), c(ddg_low, ddg_hig))
}

# CARTO Temps: teal (stabilizing) → gold → rose (destabilizing).
fill_cols <- palette_colors("Diverging.Temps")
zero_at <- (0 - fill_lim[[1]]) / (fill_lim[[2]] - fill_lim[[1]])
zero_at <- min(max(zero_at, 0.05), 0.95)
n_fill <- length(fill_cols)
half <- (n_fill - 1L) / 2
fill_values <- c(
    seq(0, zero_at, length.out = half + 1L),
    seq(zero_at, 1, length.out = half + 1L)[-1]
)
algolia <- palette_colors("Brand.Algolia")
safe <- palette_colors("Qualitative.Safe")
col_hot <- algolia[[5]]
col_ok <- algolia[[3]]
col_line <- safe[[12]]

blank_bg <- theme(
    panel.grid = element_blank(),
    plot.background = element_blank(),
    panel.background = element_blank(),
    legend.background = element_blank()
)

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p_heat <- ggplot() +
    geom_tile(
        data = heat,
        aes(x = x, y = y, fill = ddg),
        width = 1,
        height = 1
    ) +
    geom_hline(
        yintercept = y_zero,
        linetype = "dashed",
        colour = algolia[[10]],
        linewidth = 0.7
    ) +
    geom_line(
        data = avg,
        aes(x = x, y = y_line),
        colour = col_line,
        linetype = "dashed",
        linewidth = 0.85
    ) +
    geom_point(
        data = avg,
        aes(x = x, y = y_line, colour = hot),
        size = 3.2
    ) +
    scale_fill_gradientn(
        colours = fill_cols,
        values = fill_values,
        limits = fill_lim,
        oob = scales::squish,
        guide = "none"
    ) +
    scale_colour_manual(
        values = c("TRUE" = col_hot, "FALSE" = col_ok),
        guide = "none"
    ) +
    scale_x_continuous(
        breaks = seq_len(n_site),
        labels = levels(df[[site_col]]),
        limits = c(0.5, n_site + 0.5),
        expand = c(0, 0)
    ) +
    scale_y_continuous(
        breaks = seq_len(n_aa),
        labels = aa_cols,
        limits = c(0.5, n_aa + 0.5),
        expand = c(0, 0),
        sec.axis = sec_axis(
            transform = to_ddg,
            name = config$labels$y_right,
            breaks = pretty(c(ddg_low, ddg_hig), n = 6)
        )
    ) +
    labs(x = config$labels$x, y = config$labels$y) +
    theme_prism() +
    blank_bg +
    theme(
        plot.margin = margin(2, 4, 4, 4),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    )

# Same x panel as the heatmap so aplot xlim2 / patchwork keep the
# colorbar the same length as the tile grid. Tick labels are ΔΔG.
x_lim <- c(0.5, n_site + 0.5)
fill_breaks <- pretty(fill_lim, n = 6)
fill_breaks <- fill_breaks[
    fill_breaks >= fill_lim[[1]] & fill_breaks <= fill_lim[[2]]
]
n_bar <- 256
p_bar <- ggplot(
    data.frame(
        x = seq(x_lim[[1]], x_lim[[2]], length.out = n_bar),
        y = 1,
        fill = seq(fill_lim[[1]], fill_lim[[2]], length.out = n_bar)
    ),
    aes(x = x, y = y, fill = fill)
) +
    geom_raster(interpolate = TRUE) +
    scale_fill_gradientn(
        colours = fill_cols,
        values = fill_values,
        limits = fill_lim,
        guide = "none"
    ) +
    scale_x_continuous(
        name = config$labels$fill,
        breaks = linmap(fill_breaks, fill_lim, x_lim),
        labels = fill_breaks,
        limits = x_lim,
        expand = c(0, 0),
        position = "top"
    ) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_prism() +
    blank_bg +
    theme(
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x.bottom = element_blank()
    )

p_gap <- ggplot() +
    theme_void() +
    theme(
        plot.background = element_blank(),
        panel.background = element_blank()
    )

p <- as.patchwork(
    insert_top(insert_top(p_heat, p_gap, height = 0.05), p_bar, height = 0.11)
) +
    plot_annotation(
        theme = theme(plot.background = element_rect(fill = NA, colour = NA))
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
suppressWarnings(save_ggplot(p, io$output, width = 9.2, height = 5.8))

