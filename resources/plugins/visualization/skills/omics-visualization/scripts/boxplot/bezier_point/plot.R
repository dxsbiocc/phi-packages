#!/usr/bin/env Rscript

# Template-ID: boxplot-bezier-point
#
# Purpose:
#   Draw faceted paired boxplots with Bezier curves connecting
#   observations that share a pair_id within each facet.
#
# Inputs:
#   A table with one row per unit-by-condition. Default example:
#     - group: paired group (two levels)
#     - exprs: numeric value
#     - gene: facet variable
#     - pair_id: pairing id shared by matched observations
#
# Output:
#   A PDF, PNG, or SVG faceted boxplot with Bezier links.
#
# Dependencies:
#   ggplot2, readr, ggforce, ggpubr
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (column names and labels).
#   Edit DATA PREPARATION to change facet order or validation.
#   Edit PLOT only when the chart geometry or styling must change.
#
# Scientific assumptions:
#   Rows that share the same pair_id within a facet are matched
#   observations; pairing is not inferred from row order.
#   Bezier curves are a display transform, not a statistical model.
#   Significance labels compare the two groups within each facet.

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
        x = "group",
        y = "exprs",
        facet = "gene",
        pair_id = "pair_id"
    ),
    labels = list(
        title = "",
        x = "",
        y = "Gene expression level"
    )
)

load_packages(c("ggplot2", "readr", "ggforce", "ggpubr"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- as.data.frame(read_table_auto(io$input), stringsAsFactors = FALSE)
require_columns(df, config$columns)
x_col <- config$columns$x
y_col <- config$columns$y
facet_col <- config$columns$facet
pair_col <- config$columns$pair_id

df[[x_col]] <- factor(df[[x_col]], levels = unique(as.character(df[[x_col]])))
df[[facet_col]] <- factor(
    as.character(df[[facet_col]]),
    levels = sort(unique(as.character(df[[facet_col]])))
)
df[[pair_col]] <- as.character(df[[pair_col]])
df[[y_col]] <- as.numeric(df[[y_col]])
if (any(!is.finite(df[[y_col]]))) {
    stop("y values must be numeric and finite.", call. = FALSE)
}
if (anyNA(df[[pair_col]]) || any(!nzchar(df[[pair_col]]))) {
    stop("pair_id must be non-empty for every row.", call. = FALSE)
}
df$x_plot <- as.numeric(df[[x_col]])
facet_levels <- levels(df[[facet_col]])
x_lvls <- levels(df[[x_col]])
if (length(x_lvls) < 2L) {
    stop("Need at least two x / group levels to draw pair links.", call. = FALSE)
}

# Bezier control points pairing rows by pair_id within each facet.
bezier_rows <- list()
k <- 1L
for (facet in facet_levels) {
    subf <- df[df[[facet_col]] == facet, , drop = FALSE]
    pair_ids <- unique(subf[[pair_col]])
    for (pid in pair_ids) {
        sub <- subf[subf[[pair_col]] == pid, , drop = FALSE]
        if (anyDuplicated(as.character(sub[[x_col]]))) {
            stop(
                sprintf(
                    "pair_id '%s' in facet '%s' has duplicate group levels.",
                    pid, facet
                ),
                call. = FALSE
            )
        }
        sub <- sub[order(sub$x_plot), , drop = FALSE]
        if (nrow(sub) < 2L) {
            next
        }
        for (i in seq_len(nrow(sub) - 1L)) {
            x1 <- sub$x_plot[[i]]
            x2 <- sub$x_plot[[i + 1L]]
            y1 <- sub[[y_col]][[i]]
            y2 <- sub[[y_col]][[i + 1L]]
            bezier_rows[[k]] <- data.frame(
                bx = c(x1, x1 + 0.3, x2 - 0.3, x2),
                by = c(y1, y1, y2, y2),
                bindex = paste(facet, pid, i, sep = "_"),
                stringsAsFactors = FALSE
            )
            bezier_rows[[k]][[facet_col]] <- facet
            k <- k + 1L
        }
    }
}
if (!length(bezier_rows)) {
    stop(
        "No complete pairs found. Each pair_id within a facet needs >= 2 group levels.",
        call. = FALSE
    )
}
bezier_data <- do.call(rbind, bezier_rows)
bezier_data[[facet_col]] <- factor(
    bezier_data[[facet_col]],
    levels = facet_levels
)

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
colour_values <- c("#ed6ca4", "#fbb05b", "#acd372", "#7bc4e2")

p <- ggplot(df, aes(
    x = x_plot,
    y = .data[[y_col]],
    colour = .data[[x_col]]
)) +
    geom_boxplot(outlier.shape = NA, show.legend = FALSE) +
    geom_bezier(
        data = bezier_data,
        aes(x = bx, y = by, group = bindex),
        colour = "#caccd1",
        alpha = 0.7,
        inherit.aes = FALSE
    ) +
    geom_point(
        position = position_dodge2(width = 0.2),
        size = 3,
        show.legend = FALSE,
        alpha = 0.7
    ) +
    stat_compare_means(
        aes(group = .data[[x_col]]),
        label = "p.signif",
        hide.ns = TRUE,
        show.legend = FALSE,
        label.x = 1.5
    ) +
    scale_colour_manual(values = expand_palette(colour_values, length(x_lvls))) +
    scale_x_continuous(
        breaks = seq_along(x_lvls),
        labels = x_lvls
    ) +
    facet_wrap(
        vars(.data[[facet_col]]),
        nrow = 1,
        scales = "free_x"
    ) +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y
    ) +
    theme_bw() +
    theme(
        plot.background = element_blank(),
        axis.text = element_text(face = "bold", size = 8),
        axis.title.y = element_text(size = 10, face = "bold"),
        strip.text = element_text(face = "bold"),
        strip.background = element_rect(fill = "#a0c44387"),
        panel.grid = element_blank(),
        panel.background = element_blank()
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output)
