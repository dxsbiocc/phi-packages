#!/usr/bin/env Rscript

# Template-ID: scatter-volcano-inset
#
# Purpose:
#   Draw a volcano (effect vs -log10 p) with a ggmagnify callout:
#   a dashed target window, projection lines, and a labelled inset.
#
# Inputs:
#   One row per tested feature. Default example reuses the
#   scatter-volcano table (ceiling q-values dropped):
#     - symbol: label
#     - log2FC: signed effect
#     - qvalue: supplied adjusted p-value
#
# Output:
#   A PDF, PNG, or SVG volcano plot with a zoom inset.
#
# Dependencies:
#   ggplot2, readr, ggprism, ggrepel, ggmagnify
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (columns, cutoffs, from/to,
#   main_labels). Edit DATA PREPARATION to change how significance
#   is split. Edit PLOT only when the chart geometry must change.
#   The callout is geom_magnify(from, to, plot = p_zoom). Do not
#   replace it with annotation_custom() or a second template id.
#
# Scientific assumptions:
#   Effect and p are supplied results. This script does not run a
#   differential test. Colour is a display split: p versus p_cutoff
#   and |x| versus fc_line (Down / Up / Non-significant).
#   Vertical dashed lines are at ±fc_line, not at x = 0.
#   from/to are CONFIG windows, not a peak caller.
#   -log10(p) is a display transform; p is floored at 1e-30.

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
        x = "log2FC",
        p = "qvalue",
        label = "symbol"
    ),
    p_cutoff = 0.05,
    fc_line = 1,
    xlim = c(-8, 10),
    ylim = c(0, 24),
    # geom_magnify: from = target window, to = inset box (data units).
    from = c(-2.25, -1.05, 1.30, 6.30),
    to = c(3.2, 8.0, 8.2, 23.2),
    main_labels = c("Sult5a1", "Cyp17a1", "Tk1", "Chil1"),
    labels = list(
        title = "",
        x = "log2FC",
        y = expression(-log[10](italic(q))),
        colour = NULL
    )
)

load_packages(c("ggplot2", "readr", "ggprism", "ggrepel", "ggmagnify"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
require_columns(df, config$columns)

x_col <- config$columns$x
p_col <- config$columns$p
lab_col <- config$columns$label
p_cutoff <- config$p_cutoff
fc_line <- abs(as.numeric(config$fc_line)[[1]])
y_line <- -log10(p_cutoff)
from <- as.numeric(config$from)
to <- as.numeric(config$to)
main_ids <- as.character(config$main_labels)

df[[x_col]] <- as.numeric(df[[x_col]])
df[[p_col]] <- pmax(as.numeric(df[[p_col]]), 1e-30)
ok <- is.finite(df[[x_col]]) & is.finite(df[[p_col]]) & df[[p_col]] > 0
if (any(!ok)) {
    message("Dropped ", sum(!ok), " row(s) with invalid effect or p.")
    df <- df[ok, , drop = FALSE]
}
if (!nrow(df)) {
    stop("No rows left to plot.", call. = FALSE)
}

df$neg_log10 <- -log10(df[[p_col]])
lab_down <- "Down"
lab_up <- "Up"
lab_ns <- "Non-significant"
df$sig <- lab_ns
df$sig[df[[p_col]] <= p_cutoff & df[[x_col]] <= -fc_line] <- lab_down
df$sig[df[[p_col]] <= p_cutoff & df[[x_col]] >= fc_line] <- lab_up
df$sig <- factor(df$sig, levels = c(lab_down, lab_ns, lab_up))

in_from <- df[[x_col]] >= min(from[1:2]) & df[[x_col]] <= max(from[1:2]) &
    df$neg_log10 >= min(from[3:4]) & df$neg_log10 <= max(from[3:4])
main_lab <- df[df[[lab_col]] %in% main_ids, , drop = FALSE]
zoom_lab <- df[
    df$sig != lab_ns & in_from & !df[[lab_col]] %in% main_ids,
    ,
    drop = FALSE
]
ns_df <- df[df$sig == lab_ns, , drop = FALSE]
down_df <- df[df$sig == lab_down, , drop = FALSE]
up_df <- df[df$sig == lab_up, , drop = FALSE]

safe <- palette_colors("Qualitative.Safe")
col_down <- safe[[4]]
col_up <- safe[[2]]
col_ns <- safe[[12]]
col_zoom <- palette_colors("Qualitative.Bold")[[5]]
fill_vals <- c(col_down, col_ns, col_up)
names(fill_vals) <- levels(df$sig)
colour_name <- config$labels$colour

shared_theme <- theme_prism() +
    theme(
        panel.grid = element_blank(),
        plot.background = element_blank(),
        panel.background = element_blank(),
        legend.background = element_blank()
    )

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p_zoom <- ggplot(df, aes(x = .data[[x_col]], y = neg_log10)) +
    geom_hline(
        yintercept = y_line,
        linetype = "dashed",
        colour = col_ns,
        linewidth = 0.45
    ) +
    geom_vline(
        xintercept = c(-fc_line, fc_line),
        linetype = "dashed",
        colour = col_ns,
        linewidth = 0.45
    ) +
    geom_point(
        data = ns_df,
        aes(colour = sig),
        size = 1.6,
        alpha = 0.35,
        stroke = 0
    ) +
    geom_point(
        data = down_df,
        aes(colour = sig),
        size = 2.6,
        alpha = 0.8,
        stroke = 0
    ) +
    geom_point(
        data = up_df,
        aes(colour = sig),
        size = 2.6,
        alpha = 0.8,
        stroke = 0
    ) +
    scale_colour_manual(values = fill_vals, drop = FALSE, guide = "none") +
    geom_text_repel(
        data = zoom_lab,
        aes(label = .data[[lab_col]], colour = sig),
        size = 3,
        seed = 1,
        box.padding = 0.3,
        point.padding = 0.3,
        segment.colour = col_ns,
        segment.size = 0.3,
        max.overlaps = Inf,
        show.legend = FALSE
    ) +
    scale_x_continuous(limits = config$xlim, expand = c(0, 0)) +
    scale_y_continuous(limits = config$ylim, expand = c(0, 0)) +
    shared_theme +
    theme(
        axis.title = element_blank(),
        legend.position = "none"
    )

p <- ggplot(df, aes(x = .data[[x_col]], y = neg_log10)) +
    geom_hline(
        yintercept = y_line,
        linetype = "dashed",
        colour = col_ns,
        linewidth = 0.45
    ) +
    geom_vline(
        xintercept = c(-fc_line, fc_line),
        linetype = "dashed",
        colour = col_ns,
        linewidth = 0.45
    ) +
    geom_point(
        data = ns_df,
        aes(colour = sig),
        size = 1.2,
        alpha = 0.35,
        stroke = 0
    ) +
    geom_point(
        data = down_df,
        aes(colour = sig),
        size = 2.2,
        alpha = 0.75,
        stroke = 0
    ) +
    geom_point(
        data = up_df,
        aes(colour = sig),
        size = 2.2,
        alpha = 0.75,
        stroke = 0
    ) +
    scale_colour_manual(values = fill_vals, drop = FALSE, name = colour_name) +
    guides(colour = guide_legend(override.aes = list(alpha = 1, size = c(4, 3, 4)))) +
    geom_text_repel(
        data = main_lab,
        aes(label = .data[[lab_col]], colour = sig),
        size = 4,
        fontface = "italic",
        seed = 1,
        box.padding = 0.3,
        point.padding = 0.3,
        segment.colour = col_ns,
        segment.size = 0.3,
        max.overlaps = Inf,
        show.legend = FALSE
    ) +
    scale_x_continuous(limits = config$xlim, expand = c(0, 0)) +
    scale_y_continuous(limits = config$ylim, expand = c(0, 0)) +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y
    ) +
    coord_cartesian(clip = "off") +
    theme_prism() +
    theme(
        panel.grid = element_blank(),
        plot.background = element_blank(),
        panel.background = element_blank(),
        legend.background = element_blank(),
        plot.margin = margin(t = 0.4, r = 0.8, b = 0.5, l = 0.5, unit = "cm"),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.justification = "center"
    ) +
    geom_magnify(
        from = from,
        to = to,
        plot = p_zoom,
        colour = col_zoom,
        target.linetype = 2,
        proj.linetype = 2,
        axes = "xy",
        expand = 0
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output, width = 9.2, height = 6.4)
