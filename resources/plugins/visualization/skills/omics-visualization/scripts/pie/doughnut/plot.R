#!/usr/bin/env Rscript

# Template-ID: pie-doughnut
#
# Purpose:
#   Draw a doughnut (ring) chart of part-to-whole composition.
#
# Inputs:
#   A two-column table (tsv/csv). Default example:
#     - name: slice label
#     - value: numeric size
#
# Output:
#   A PDF, PNG, or SVG doughnut chart.
#
# Dependencies:
#   ggplot2, readr, ggprism
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (column names and labels).
#   Edit DATA PREPARATION to change category order or reshape the table.
#   Edit PLOT only when the chart geometry or styling must change.
#
# Scientific assumptions:
#   Each row is one mutually exclusive part of a whole.
#   The inner hole is a display offset, not missing data.
#   Values are treated as already summarized, not as raw replicates.
#
# Segment shapes:
#   segment_shape = "standard" draws the default polar ring with a legend.
#   segment_shape = "chevron" draws interlocking arrow-shaped segments
#   pointing clockwise, labels written along the ring and an optional
#   centre title (example_chevron.tsv). Each segment's tip and notch have
#   the same depth, so angular span still equals value / sum. Equal values
#   give a conceptual cycle diagram.

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
        x = "name",
        y = "value"
    ),
    labels = list(
        title = "Doughnut Chart",
        x = "",
        y = ""
    ),
    # "standard" or "chevron"
    segment_shape = "standard",
    # NULL keeps the default colours; otherwise a catalog palette id
    palette = NULL,
    # The settings below apply to segment_shape = "chevron"
    hole = 0.52,                # inner radius / outer radius
    chevron_depth_deg = 9,      # how far each tip pushes into the next segment
    gap_deg = 1.6,              # white gap between segments
    start_degree = 0,           # first segment starts here, clockwise from 12 o'clock
    label_wrap = 14,            # wrap labels at this many characters
    label_size = 5.4,
    flip_bottom_labels = TRUE,  # keep lower-half labels upright
    value_labels = FALSE,       # append the percentage under each label
    center_label = NULL,        # text in the hole; use "\n" for line breaks
    center_size = 8
)

load_packages(c("ggplot2", "readr", "ggprism"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
require_columns(df, config$columns)
x_col <- config$columns$x
y_col <- config$columns$y
df[[x_col]] <- factor(df[[x_col]], levels = unique(df[[x_col]]))
df[[y_col]] <- as.numeric(df[[y_col]])

fill_values <- if (is.null(config$palette)) {
    expand_palette(
        c("#f29556",
        "#84afa2",
        "#f24c7c",
        "#29c4d0",
        "#72c427"),
        nlevels(df[[x_col]])
    )
} else {
    palette_colors(config$palette, nlevels(df[[x_col]]))
}

is_chevron <- identical(config$segment_shape, "chevron")
if (is_chevron) {
    total <- sum(df[[y_col]])
    span <- 360 * df[[y_col]] / total
    t_start <- config$start_degree + cumsum(span) - span + config$gap_deg / 2
    t_end <- t_start + span - config$gap_deg
    depth <- min(config$chevron_depth_deg, 0.45 * min(t_end - t_start))
    r_in <- config$hole
    r_mid <- (r_in + 1) / 2
    polar_xy <- function(r, deg) {
        data.frame(x = r * sin(deg * pi / 180), y = r * cos(deg * pi / 180))
    }

    shapes <- do.call(rbind, lapply(seq_len(nrow(df)), function(i) {
        arc <- seq(t_start[[i]], t_end[[i]], length.out = 60)
        pts <- rbind(
            polar_xy(1, arc),
            polar_xy(r_mid, t_end[[i]] + depth),
            polar_xy(r_in, rev(arc)),
            polar_xy(r_mid, t_start[[i]] + depth)
        )
        pts$segment <- df[[x_col]][[i]]
        pts
    }))

    label_deg <- (t_start + t_end) / 2 + depth / 2
    label_text <- vapply(as.character(df[[x_col]]), function(s) {
        paste(strwrap(s, width = config$label_wrap), collapse = "\n")
    }, character(1), USE.NAMES = FALSE)
    if (isTRUE(config$value_labels)) {
        label_text <- paste0(label_text, "\n",
            sprintf("%.1f%%", 100 * df[[y_col]] / total))
    }
    text_angle <- -label_deg
    if (isTRUE(config$flip_bottom_labels)) {
        lower <- (label_deg %% 360) > 90 & (label_deg %% 360) < 270
        text_angle[lower] <- text_angle[lower] + 180
    }
    labels_df <- cbind(polar_xy(r_mid, label_deg),
        label = label_text, angle = text_angle)
}

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
blank_theme <- theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5),
    plot.background = element_blank(),
    legend.background = element_blank()
)

if (is_chevron) {
    p <- ggplot() +
        geom_polygon(
            data = shapes,
            aes(x = x, y = y, group = segment, fill = segment),
            colour = NA
        ) +
        geom_text(
            data = labels_df,
            aes(x = x, y = y, label = label, angle = angle),
            size = config$label_size, lineheight = 0.9, colour = "#222222"
        ) +
        scale_fill_manual(values = fill_values, guide = "none") +
        coord_equal(xlim = c(-1.02, 1.02), ylim = c(-1.02, 1.02)) +
        labs(title = config$labels$title, x = NULL, y = NULL) +
        theme_prism() +
        blank_theme
    if (!is.null(config$center_label)) {
        p <- p + annotate("text", x = 0, y = 0, label = config$center_label,
            size = config$center_size, lineheight = 1, colour = "#222222")
    }
} else {
    p <- ggplot(df, aes(x = 2, y = .data[[y_col]], fill = .data[[x_col]])) +
        geom_col(width = 1, colour = "white") +
        scale_fill_manual(values = fill_values) +
        coord_polar(theta = "y") +
        xlim(0.5, 2.5) +
        labs(title = config$labels$title, x = NULL, y = NULL, fill = "") +
        theme_prism() +
        blank_theme
}

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output, width = 8, height = 8)
