#!/usr/bin/env Rscript

# Template-ID: line-group
#
# Purpose:
#   Draw a grouped line chart of a numeric value across ordered
#   categories, with a line and points per group.
#
# Inputs:
#   A table with one row per category-by-group. Default example:
#     - day: ordered category
#     - value: numeric value
#     - category: group
#   style = "growth" (example_growth.tsv) accepts either replicate rows
#   (one row per subject per time point; mean and SD / SEM are
#   computed as display summaries) or one summarised row per
#   time-by-group with columns$error holding the supplied error.
#   Optional sidecar brackets.tsv: group1, group2, x, label. The label
#   (e.g. "**") is supplied; no test is run here.
#
# Output:
#   A PDF, PNG, or SVG grouped line chart.
#
# Dependencies:
#   ggplot2, readr
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (column names and labels).
#   Edit DATA PREPARATION to change category or group order.
#   Edit PLOT only when the chart geometry or styling must change.
#
# Scientific assumptions:
#   Each row is one sequential observation in a group.
#   Category order in the file is treated as the x-axis order.
#   Connecting points with a line assumes an ordered relationship.

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
    # "basic": one line per group from one value per point.
    # "growth": mean +/- error per point (e.g. tumour volume over days),
    #           in-panel legend, rotated x labels, supplied brackets.
    style = "basic",
    columns = list(
        x = "day",
        y = "value",
        colour = "category",
        error = NULL  # growth only: supplied error column (y is then the mean)
    ),
    labels = list(
        title = "Lines by Category",
        x = "Country",
        y = "Income",
        y_expr = NULL  # plotmath string overriding y, e.g. "Tumour~volume~(mm^3)"
    ),
    growth = list(
        error = "sem",  # "sem" | "sd" | "none" (replicate input only)
        palette = c(
            "Qualitative.Safe:12", "Qualitative.Bold:8",
            "Qualitative.Bold:1", "Qualitative.Bold:5"
        ),
        point_size = 4.2,
        point_alpha = 0.8,
        line_width = 0.5,
        errorbar_width = 0.18,
        errorbar_linewidth = 0.5,
        legend_inside = c(0.03, 0.98),  # NULL puts the legend on the right
        x_angle = 45,
        brackets_file = "brackets.tsv",
        bracket_offset = 0.32,  # x distance right of the point, in category units
        bracket_tick = 0.08,
        label_size = 6,
        base_size = 14,
        width = 5.5,
        height = 5
    )
)

load_packages(c("ggplot2", "readr"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
require_columns(df, config$columns)
x_col <- config$columns$x
y_col <- config$columns$y
colour_col <- config$columns$colour
df[[x_col]] <- factor(df[[x_col]], levels = unique(df[[x_col]]))
df[[colour_col]] <- factor(df[[colour_col]], levels = unique(df[[colour_col]]))
df[[y_col]] <- as.numeric(df[[y_col]])

pal_ref <- function(ref) {
    parts <- strsplit(ref, ":", fixed = TRUE)[[1]]
    cols <- palette_colors(parts[[1]])
    if (length(parts) > 1L) cols[[as.integer(parts[[2]])]] else cols
}

sidecar_path <- function(name) {
    if (is.null(name)) return(NULL)
    path <- if (file.exists(name)) name else file.path(dirname(io$input), name)
    if (file.exists(path)) path else NULL
}

y_title <- if (is.null(config$labels$y_expr)) {
    config$labels$y
} else {
    parse(text = config$labels$y_expr)[[1]]
}

if (identical(config$style, "growth")) {
    gc <- config$growth
    err_col <- config$columns$error
    if (is.null(err_col)) {
        if (!gc$error %in% c("sem", "sd", "none")) {
            stop("growth$error must be sem, sd or none", call. = FALSE)
        }
        keys <- unique(df[c(x_col, colour_col)])
        sm <- do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
            v <- df[[y_col]][df[[x_col]] == keys[[x_col]][i] &
                df[[colour_col]] == keys[[colour_col]][i]]
            v <- v[is.finite(v)]
            s <- if (length(v) > 1L) stats::sd(v) else 0
            data.frame(
                x = keys[[x_col]][i], g = keys[[colour_col]][i], mean = mean(v),
                err = switch(gc$error, sd = s, sem = s / sqrt(length(v)), none = 0)
            )
        }))
    } else {
        sm <- data.frame(
            x = df[[x_col]], g = df[[colour_col]],
            mean = df[[y_col]], err = as.numeric(df[[err_col]])
        )
    }
    sm$x <- factor(sm$x, levels = levels(df[[x_col]]))
    sm$g <- factor(sm$g, levels = levels(df[[colour_col]]))
    sm <- sm[order(sm$g, sm$x), ]
    sm$xi <- as.integer(sm$x)

    groups <- levels(sm$g)
    cols <- unlist(lapply(gc$palette, pal_ref), use.names = FALSE)
    cols <- setNames(expand_palette(cols, length(groups)), groups)

    br <- NULL
    br_path <- sidecar_path(gc$brackets_file)
    if (!is.null(br_path)) {
        b <- read_table_auto(br_path)
        require_columns(b, list("group1", "group2", "x", "label"))
        br <- do.call(rbind, lapply(seq_len(nrow(b)), function(i) {
            y1 <- sm$mean[sm$g == b$group1[i] & sm$x == b$x[i]]
            y2 <- sm$mean[sm$g == b$group2[i] & sm$x == b$x[i]]
            if (!length(y1) || !length(y2)) {
                stop(sprintf(
                    "Bracket %s vs %s at %s has no matching point.",
                    b$group1[i], b$group2[i], b$x[i]
                ), call. = FALSE)
            }
            xb <- match(b$x[i], levels(sm$x)) + gc$bracket_offset
            data.frame(x = xb, y1 = y1, y2 = y2, label = b$label[i])
        }))
    }

    x_hi <- length(levels(sm$x)) + if (is.null(br)) 0.4 else gc$bracket_offset + 0.45
    p <- ggplot(sm, aes(x = xi, y = mean, colour = g, group = g)) +
        geom_line(linewidth = gc$line_width) +
        geom_errorbar(
            aes(ymin = mean - err, ymax = mean + err),
            width = gc$errorbar_width, linewidth = gc$errorbar_linewidth
        ) +
        geom_point(size = gc$point_size, alpha = gc$point_alpha, stroke = 0) +
        scale_colour_manual(values = cols, name = NULL) +
        scale_x_continuous(
            breaks = seq_along(levels(sm$x)), labels = levels(sm$x),
            limits = c(0.6, x_hi), expand = c(0, 0)
        ) +
        scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
        labs(title = config$labels$title, x = config$labels$x, y = y_title) +
        theme_classic(base_size = gc$base_size) +
        theme(
            plot.title = element_text(hjust = 0.5, face = "bold"),
            axis.text = element_text(colour = "black"),
            axis.text.x = element_text(
                angle = gc$x_angle,
                hjust = if (gc$x_angle) 1 else 0.5,
                vjust = if (gc$x_angle) 1 else 0.5
            ),
            axis.ticks = element_line(colour = "black"),
            legend.background = element_blank(),
            legend.key = element_blank(),
            plot.background = element_blank(),
            panel.background = element_blank()
        )
    if (is.null(gc$legend_inside)) {
        p <- p + theme(legend.position = "right")
    } else {
        p <- p + theme(
            legend.position = "inside",
            legend.position.inside = gc$legend_inside,
            legend.justification = c(0, 1)
        )
    }
    if (!is.null(br)) {
        tk <- gc$bracket_tick
        p <- p +
            geom_segment(
                data = br, aes(x = x, xend = x, y = y1, yend = y2),
                inherit.aes = FALSE, linewidth = 0.6
            ) +
            geom_segment(
                data = rbind(
                    data.frame(x = br$x, y = br$y1),
                    data.frame(x = br$x, y = br$y2)
                ),
                aes(x = x - tk, xend = x, y = y, yend = y),
                inherit.aes = FALSE, linewidth = 0.6
            ) +
            geom_text(
                data = br, aes(x = x + 0.04, y = (y1 + y2) / 2, label = label),
                inherit.aes = FALSE, hjust = 0, vjust = 0.75,
                size = gc$label_size, fontface = "bold"
            )
    }
    if (is.null(config$labels$title)) p <- p + theme(plot.title = element_blank())
    save_ggplot(p, io$output, width = gc$width, height = gc$height)
    quit(save = "no")
}

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p <- ggplot(df, aes(
    x = .data[[x_col]],
    y = .data[[y_col]],
    colour = .data[[colour_col]]
)) +
    geom_line(aes(group = .data[[colour_col]])) +
    geom_point(colour = "white", size = 3) +
    geom_point(size = 2) +
    scale_colour_viridis_d() +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y
    ) +
    theme_minimal() +
    theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "bottom",
        plot.background = element_blank()
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output)
