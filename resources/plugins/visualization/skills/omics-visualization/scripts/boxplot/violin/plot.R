#!/usr/bin/env Rscript

# Template-ID: boxplot-violin
#
# Purpose:
#   Draw a violin plot of a numeric value across categories, with a
#   nested boxplot in each violin. style = "jitter" adds coloured
#   violins, a thin box, and jittered points, with an optional
#   supplied n / P annotation.
#
# Inputs:
#   A two-column table (tsv/csv). Default example:
#     - gene: category label
#     - expression: numeric value
#   Jitter example (example_jitter.tsv):
#     - group, value  (e.g. year|scenario habitat change)
#
# Output:
#   A PDF, PNG, or SVG violin plot.
#
# Dependencies:
#   ggplot2, readr, ggprism
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (column names and labels).
#   Edit DATA PREPARATION to change category order.
#   Edit PLOT only when the chart geometry or styling must change.
#   Do not run significance tests here; annotation$label is supplied.
#   Half-violin rainclouds stay boxplot-raincloud. Mirrored score
#   pairs stay boxplot-violin-butterfly.
#
# Scientific assumptions:
#   Each row is one observation of a numeric value in a category.
#   Violins estimate the density of the raw values.
#   Nested boxes summarize the same distribution.

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
    # "basic": nested box in violin (default preview, theme_prism).
    # "jitter": filled violin + thin box + jitter points, white panel.
    style = "basic",
    columns = list(
        x = "gene",
        y = "expression"
    ),
    labels = list(
        title = "",
        x = "",
        y = "Gene expression level"
    ),
    group_order = NULL,
    palettes = list(
        group = "Qualitative.Bold"
    ),
    jitter = list(
        violin_alpha = 0.35,
        violin_linewidth = 0.4,
        box_width = 0.12,
        box_linewidth = 0.45,
        point_size = 1.8,
        point_alpha = 0.75,
        jitter_width = 0.08,
        x_angle = 90,
        # Optional corner note; text is supplied (e.g. "n = 14\nP < 0.001").
        annotation = list(
            label = NULL,
            x = -Inf,
            y = Inf,
            hjust = -0.1,
            vjust = 1.2,
            size = 3.6
        ),
        base_size = 13,
        width = 7.2,
        height = 5.2
    )
)

load_packages(c("ggplot2", "readr", "ggprism"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- as.data.frame(read_table_auto(io$input), stringsAsFactors = FALSE)
require_columns(df, config$columns)
x_col <- config$columns$x
y_col <- config$columns$y
df[[y_col]] <- as.numeric(df[[y_col]])
lv <- unique(as.character(df[[x_col]]))
if (!is.null(config$group_order)) {
    lv <- c(
        intersect(as.character(config$group_order), lv),
        setdiff(lv, as.character(config$group_order))
    )
}
df[[x_col]] <- factor(df[[x_col]], levels = lv)
df <- df[is.finite(df[[y_col]]) & !is.na(df[[x_col]]), , drop = FALSE]
if (!nrow(df)) stop("No finite rows to plot.", call. = FALSE)

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
if (identical(config$style, "jitter")) {
    jc <- config$jitter
    cols <- expand_palette(
        palette_colors(config$palettes$group), nlevels(df[[x_col]])
    )
    names(cols) <- levels(df[[x_col]])
    p <- ggplot(df, aes(
        x = .data[[x_col]],
        y = .data[[y_col]],
        colour = .data[[x_col]],
        fill = .data[[x_col]]
    )) +
        geom_violin(
            alpha = jc$violin_alpha,
            linewidth = jc$violin_linewidth,
            show.legend = FALSE
        ) +
        geom_boxplot(
            width = jc$box_width,
            outlier.shape = NA,
            fill = NA,
            linewidth = jc$box_linewidth,
            show.legend = FALSE
        ) +
        geom_jitter(
            width = jc$jitter_width,
            height = 0,
            size = jc$point_size,
            alpha = jc$point_alpha,
            show.legend = FALSE
        ) +
        scale_fill_manual(values = cols, guide = "none") +
        scale_colour_manual(values = cols, guide = "none") +
        labs(
            title = if (is.null(config$labels$title) ||
                !nzchar(as.character(config$labels$title))) {
                NULL
            } else {
                config$labels$title
            },
            x = config$labels$x,
            y = config$labels$y
        ) +
        theme_classic(base_size = jc$base_size) +
        theme(
            axis.text = element_text(colour = "black"),
            axis.text.x = element_text(
                angle = jc$x_angle,
                hjust = if (jc$x_angle) 1 else 0.5,
                vjust = if (jc$x_angle) 0.5 else 0.5
            ),
            axis.ticks = element_line(colour = "black"),
            plot.background = element_rect(fill = "white", colour = NA),
            panel.background = element_rect(fill = "white", colour = NA)
        )
    ann <- jc$annotation
    if (!is.null(ann$label) && nzchar(as.character(ann$label))) {
        p <- p + annotate(
            "text",
            x = ann$x, y = ann$y, label = ann$label,
            hjust = ann$hjust, vjust = ann$vjust,
            size = ann$size, colour = "grey15"
        )
    }
    if (is.null(config$labels$title) ||
        !nzchar(as.character(config$labels$title))) {
        p <- p + theme(plot.title = element_blank())
    }
    save_ggplot(p, io$output, width = jc$width, height = jc$height)
    quit(save = "no")
}

fill_values <- c(
    "#009392", "#39b185", "#9ccb86", "#e9e29c",
    "#eeb479", "#e88471", "#cf597e"
)

p <- ggplot(df, aes(
    x = .data[[x_col]],
    y = .data[[y_col]]
)) +
    geom_violin() +
    geom_boxplot(
        aes(fill = .data[[x_col]]),
        width = 0.2,
        outlier.shape = NA,
        show.legend = FALSE
    ) +
    scale_fill_manual(name = NULL, values = fill_values) +
    guides(x = guide_prism_bracket()) +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y
    ) +
    theme_prism() +
    theme(plot.background = element_blank())

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output)
