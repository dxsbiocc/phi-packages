#!/usr/bin/env Rscript

# Template-ID: bar-communication
#
# Purpose:
#   Draw tissue-wise mirrored interaction bars (outgoing left, incoming
#   right) beside a log2 ratio heatmap. Rows are cell types and stay
#   aligned across panels.
#
# Inputs:
#   One table, one row per cell × tissue × direction:
#     cell, tissue, direction (outgoing|incoming), strength
#   Cell order is first-appearance order (first cell at the top).
#   Tissue order is first-appearance order. The focal tissue is the
#   last tissue unless CONFIG$focal names another level.
#
# Output:
#   A PDF, PNG, or SVG combined bar-and-heatmap figure.
#
# Dependencies:
#   ggplot2, readr, patchwork
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   Edit only CONFIG for column names, the focal tissue, the highlighted
#   cell, and the heatmap clamp.
#   Edit PLOT only when the chart geometry or styling must change.
#
# Scientific assumptions:
#   strength is a supplied interaction summary, not a score computed here.
#   log2(focal / baseline) is a display transform. Strengths must be > 0.
#   The dashed box marks CONFIG$highlight$cell. The heatmap box marks
#   that cell at highlight$direction versus highlight$baseline.
#   This script does not run CellChat or any communication model.

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
        cell = "cell",
        tissue = "tissue",
        direction = "direction",
        strength = "strength"
    ),
    focal = "T_M",
    highlight = list(
        cell = "Pericyte",
        direction = "incoming",
        baseline = "N"
    ),
    heat_limit = 3,
    labels = list(
        x = "Interaction strength",
        fill = "Scale",
        outgoing = "Outgoing",
        incoming = "Incoming"
    ),
    # Hex copied from Qualitative.Bold and Diverging.RdYlBu.
    colours = list(
        outgoing = "#f2b701",
        incoming = "#008695",
        highlight = "#d73027",
        cell = "Qualitative.Bold",
        heat = "Diverging.RdYlBu"
    )
)

load_packages(c("ggplot2", "readr", "patchwork", "scales"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- as.data.frame(read_table_auto(io$input), stringsAsFactors = FALSE)
require_columns(df, config$columns)
cell_col <- config$columns$cell
tissue_col <- config$columns$tissue
dir_col <- config$columns$direction
strength_col <- config$columns$strength

df[[cell_col]] <- as.character(df[[cell_col]])
df[[tissue_col]] <- as.character(df[[tissue_col]])
df[[dir_col]] <- tolower(trimws(as.character(df[[dir_col]])))
df[[strength_col]] <- as.numeric(df[[strength_col]])

if (!all(df[[dir_col]] %in% c("outgoing", "incoming"))) {
    stop('direction must be "outgoing" or "incoming".', call. = FALSE)
}
if (any(!is.finite(df[[strength_col]]) | df[[strength_col]] <= 0)) {
    stop("strength must be positive and finite.", call. = FALSE)
}

cell_levels <- unique(df[[cell_col]])
tissue_levels <- unique(df[[tissue_col]])
if (length(tissue_levels) < 2L) {
    stop("Need at least two tissues.", call. = FALSE)
}
focal <- config$focal
if (!nzchar(focal) || !focal %in% tissue_levels) {
    focal <- tissue_levels[[length(tissue_levels)]]
}
baselines <- setdiff(tissue_levels, focal)
key <- paste(df[[cell_col]], df[[tissue_col]], df[[dir_col]], sep = "\r")
if (anyDuplicated(key)) {
    stop("Each cell × tissue × direction must appear once.", call. = FALSE)
}
expected <- length(cell_levels) * length(tissue_levels) * 2L
if (nrow(df) != expected) {
    stop(
        "Every cell needs both directions in every tissue.",
        call. = FALSE
    )
}

# Display order: first cell at the top (highest y).
y_of <- setNames(rev(seq_along(cell_levels)), cell_levels)
df$y <- unname(y_of[df[[cell_col]]])
df$x <- ifelse(
    df[[dir_col]] == "outgoing",
    -df[[strength_col]],
    df[[strength_col]]
)

strength_max <- max(df[[strength_col]])
lim <- max(pretty(c(0, strength_max * 1.05)))
if (lim <= strength_max) {
    lim <- lim + max(diff(pretty(c(0, strength_max))))
}
x_breaks <- pretty(c(-lim, lim))
x_breaks <- x_breaks[abs(x_breaks) <= lim]

heat_rows <- list()
k <- 1L
for (direction in c("outgoing", "incoming")) {
    for (baseline in baselines) {
        for (cell in cell_levels) {
            num <- df[[strength_col]][
                df[[cell_col]] == cell &
                    df[[tissue_col]] == focal &
                    df[[dir_col]] == direction
            ]
            den <- df[[strength_col]][
                df[[cell_col]] == cell &
                    df[[tissue_col]] == baseline &
                    df[[dir_col]] == direction
            ]
            heat_rows[[k]] <- data.frame(
                cell = cell,
                y = unname(y_of[[cell]]),
                direction = direction,
                baseline = baseline,
                contrast = paste(focal, baseline, sep = "/"),
                value = log2(num / den),
                stringsAsFactors = FALSE
            )
            k <- k + 1L
        }
    }
}
heat <- do.call(rbind, heat_rows)
if (any(!is.finite(heat$value))) {
    stop("log2 ratio is not finite. Check zero or missing strengths.", call. = FALSE)
}
contrast_levels <- unique(heat$contrast)
# Column order: all outgoing contrasts, then all incoming contrasts.
heat$x <- match(paste(heat$direction, heat$contrast), unique(paste(heat$direction, heat$contrast)))

hi_cell <- config$highlight$cell
hi_dir <- tolower(config$highlight$direction)
hi_base <- config$highlight$baseline
if (!hi_cell %in% cell_levels) {
    stop(paste("Highlight cell not in the table:", hi_cell), call. = FALSE)
}
if (!hi_dir %in% c("outgoing", "incoming")) {
    stop('highlight$direction must be "outgoing" or "incoming".', call. = FALSE)
}
if (!hi_base %in% baselines) {
    stop(paste("Highlight baseline must be a non-focal tissue:", hi_base), call. = FALSE)
}
hi_y <- unname(y_of[[hi_cell]])
hi_contrast <- paste(focal, hi_base, sep = "/")
hi_x <- heat$x[heat$cell == hi_cell & heat$direction == hi_dir & heat$baseline == hi_base][[1]]

y_limits <- c(0.45, length(cell_levels) + 0.55)
cell_pal <- palette_colors(config$colours$cell, n = length(cell_levels))
names(cell_pal) <- cell_levels
heat_pal <- rev(palette_colors(config$colours$heat))
heat_lim <- config$heat_limit

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
panel_theme <- function(show_x_title = FALSE) {
    theme_classic(base_size = 11) +
        theme(
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.line.y = element_blank(),
            axis.line.x = element_line(linewidth = 0.3, colour = "black"),
            axis.ticks.x = element_line(linewidth = 0.3),
            axis.text.x = element_text(size = 8, colour = "black"),
            axis.title.x = element_text(
                size = 11,
                colour = if (show_x_title) "black" else NA
            ),
            plot.title = element_text(size = 12, hjust = 0.5, colour = "black"),
            plot.background = element_rect(fill = "white", colour = NA),
            panel.background = element_rect(fill = "white", colour = NA),
            plot.margin = margin(8, 6, 4, 6),
            legend.position = "none"
        )
}

mirror_panel <- function(tissue, show_x_title) {
    sub <- df[df[[tissue_col]] == tissue, , drop = FALSE]
    ggplot(sub, aes(x = x, y = y, fill = .data[[dir_col]])) +
        geom_col(width = 0.68, orientation = "y") +
        geom_vline(xintercept = 0, colour = "#555555", linewidth = 0.35) +
        annotate(
            "rect",
            xmin = -lim, xmax = lim,
            ymin = hi_y - 0.5, ymax = hi_y + 0.5,
            fill = NA, colour = config$colours$highlight,
            linetype = "dashed", linewidth = 0.45
        ) +
        scale_fill_manual(
            values = c(
                outgoing = config$colours$outgoing,
                incoming = config$colours$incoming
            )
        ) +
        scale_x_continuous(
            limits = c(-lim, lim),
            breaks = x_breaks,
            labels = function(x) as.character(abs(x))
        ) +
        scale_y_continuous(limits = y_limits, expand = c(0, 0)) +
        labs(title = tissue, x = config$labels$x, y = NULL) +
        panel_theme(show_x_title)
}

name_df <- data.frame(
    cell = cell_levels,
    y = unname(y_of),
    stringsAsFactors = FALSE
)
p_names <- ggplot(name_df, aes(y = y)) +
    geom_point(aes(x = 1.15, colour = cell), size = 2.4) +
    geom_text(aes(x = 1.02, label = cell, colour = cell), hjust = 1, size = 3.1) +
    scale_colour_manual(values = cell_pal, guide = "none") +
    scale_x_continuous(limits = c(-0.15, 1.35), expand = c(0, 0)) +
    scale_y_continuous(limits = y_limits, expand = c(0, 0)) +
    labs(title = " ", x = " ") +
    panel_theme(TRUE) +
    theme(
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(colour = NA),
        plot.margin = margin(8, 0, 4, 4)
    )

bar_plots <- lapply(seq_along(tissue_levels), function(i) {
    mirror_panel(tissue_levels[[i]], show_x_title = i == ceiling(length(tissue_levels) / 2))
})

x_labels <- vapply(sort(unique(heat$x)), function(i) {
    contrast <- unique(heat$contrast[heat$x == i])
    direction <- unique(heat$direction[heat$x == i])
    first_of_group <- min(heat$x[heat$direction == direction])
    prefix <- if (i == first_of_group) {
        if (direction == "outgoing") "Out: " else "In: "
    } else {
        ""
    }
    paste0(prefix, contrast)
}, character(1))
p_heat <- ggplot(heat, aes(x = x, y = y, fill = value)) +
    geom_tile(width = 0.96, height = 0.92) +
    annotate(
        "rect",
        xmin = hi_x - 0.5, xmax = hi_x + 0.5,
        ymin = hi_y - 0.5, ymax = hi_y + 0.5,
        fill = NA, colour = config$colours$highlight, linewidth = 0.7
    ) +
    scale_fill_gradientn(
        colours = heat_pal,
        limits = c(-heat_lim, heat_lim),
        oob = scales::squish,
        name = config$labels$fill,
        breaks = c(-heat_lim, 0, heat_lim)
    ) +
    scale_x_continuous(
        breaks = seq_along(x_labels),
        labels = unname(x_labels),
        expand = c(0, 0)
    ) +
    scale_y_continuous(limits = y_limits, expand = c(0, 0)) +
    labs(title = " ", x = " ", y = NULL) +
    panel_theme(TRUE) +
    theme(
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size = 8, colour = "black", angle = 45, hjust = 1),
        axis.title.x = element_text(colour = NA),
        legend.position = "right",
        legend.title = element_text(size = 10),
        legend.key.height = unit(0.55, "in"),
        plot.margin = margin(8, 4, 4, 6)
    )

p_main <- patchwork::wrap_plots(
    c(list(p_names), bar_plots, list(p_heat)),
    nrow = 1,
    widths = c(1.45, rep(1.15, length(tissue_levels)), 0.95)
)
p_legend <- ggplot() +
    annotate(
        "rect", xmin = 0, xmax = 0.28, ymin = 0.25, ymax = 0.75,
        fill = config$colours$outgoing, colour = NA
    ) +
    annotate(
        "text", x = 0.36, y = 0.5, label = config$labels$outgoing,
        hjust = 0, size = 3.4
    ) +
    annotate(
        "rect", xmin = 1.7, xmax = 1.98, ymin = 0.25, ymax = 0.75,
        fill = config$colours$incoming, colour = NA
    ) +
    annotate(
        "text", x = 2.06, y = 0.5, label = config$labels$incoming,
        hjust = 0, size = 3.4
    ) +
    scale_x_continuous(limits = c(-0.1, 4.2), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    theme_void() +
    theme(
        plot.background = element_rect(fill = "white", colour = NA),
        plot.margin = margin(2, 0, 0, 6)
    )
p <- p_legend / p_main + patchwork::plot_layout(heights = c(0.07, 1))

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output, width = 11.2, height = 5.6)
