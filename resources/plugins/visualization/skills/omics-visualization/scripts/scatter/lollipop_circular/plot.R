#!/usr/bin/env Rscript

# Template-ID: scatter-lollipop-circular
#
# Purpose:
#   Draw a circular lollipop chart with two grouping levels: outer
#   sectors (group) and lollipops (subgroup). Point radius is the
#   supplied mean; radial whiskers are the supplied se. Subgroup
#   names sit in the point centres. Numeric mean labels are not drawn.
#
# Inputs:
#   One row per lollipop. Default example:
#     - group: outer sector (T1, T2, …)
#     - subgroup: lollipop within the sector (A, B, C)
#     - mean: numeric value (radius)
#     - se: supplied standard error
#
# Output:
#   A PDF, PNG, or SVG circular lollipop chart.
#
# Dependencies:
#   ggplot2, readr, ggprism
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (column names, labels, group_gap).
#   Edit DATA PREPARATION to change group order or the y-scale ceiling.
#   Edit PLOT only when the chart geometry or styling must change.
#
# Scientific assumptions:
#   mean and se are already summarized. This script does not average
#   replicates or compute standard errors.
#   Angle is the subgroup position within its group; radius is mean.
#   Whisker length is the supplied se, not a confidence interval
#   computed here.

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
        group = "group",
        subgroup = "subgroup",
        mean = "mean",
        se = "se"
    ),
    group_gap = 1.35,
    labels = list(
        title = "Circular Lollipop Chart",
        colour = "Group"
    )
)

load_packages(c("ggplot2", "readr", "ggprism"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
require_columns(df, config$columns)

group_col <- config$columns$group
sub_col <- config$columns$subgroup
mean_col <- config$columns$mean
se_col <- config$columns$se

df[[group_col]] <- factor(df[[group_col]], levels = unique(df[[group_col]]))
df[[sub_col]] <- as.character(df[[sub_col]])
df[[mean_col]] <- as.numeric(df[[mean_col]])
df[[se_col]] <- as.numeric(df[[se_col]])

ok <- is.finite(df[[mean_col]]) & is.finite(df[[se_col]]) & df[[se_col]] >= 0
if (any(!ok)) {
    message(
        "Dropped ", sum(!ok),
        " row(s) with non-finite mean or negative/non-finite se."
    )
    df <- df[ok, , drop = FALSE]
}
if (!nrow(df)) {
    stop("No rows left to plot.", call. = FALSE)
}

group_gap <- config$group_gap
x_cursor <- 0
df$x <- NA_real_
for (g in levels(df[[group_col]])) {
    idx <- which(df[[group_col]] == g)
    n_g <- length(idx)
    df$x[idx] <- x_cursor + seq_len(n_g)
    x_cursor <- x_cursor + n_g + group_gap
}

x_min <- min(df$x) - 0.5
x_max <- max(df$x) + 0.5 + group_gap
x_span <- x_max - x_min

bands <- do.call(rbind, lapply(levels(df[[group_col]]), function(g) {
    xx <- df$x[df[[group_col]] == g]
    data.frame(
        group = g,
        xmin = min(xx) - 0.46,
        xmax = max(xx) + 0.46,
        stringsAsFactors = FALSE
    )
}))
bands$group <- factor(bands$group, levels = levels(df[[group_col]]))
bands$mid <- (bands$xmin + bands$xmax) / 2

first_mid <- bands$mid[[1]]
polar_start <- -2 * pi * (first_mid - x_min) / x_span

visual_frac <- function(x) {
    ((x - first_mid) / x_span) %% 1
}
label_angle <- function(x) {
    raw <- -360 * visual_frac(x)
    ifelse(raw < -90, raw + 180, raw)
}

df$raw_angle <- -360 * visual_frac(df$x)
df$text_angle <- ifelse(df$raw_angle < -90, df$raw_angle + 180, df$raw_angle)
bands$text_angle <- label_angle(bands$mid)

y_data_max <- max(df[[mean_col]] + df[[se_col]])
tick_max <- ceiling((y_data_max + 0.6) / 2) * 2
if (tick_max < 2) {
    tick_max <- 2
}
y_breaks <- seq(2, tick_max, by = 2)
fan_ymax <- tick_max + 0.25
band_ymin <- fan_ymax + 0.10
band_ymax <- band_ymin + 0.50
group_label_y <- band_ymax + 0.85
tick_line_x <- x_min
tick_dat <- data.frame(
    x = x_max - 0.16,
    y = y_breaks,
    label = y_breaks
)

fill_cols <- palette_colors(
    "Qualitative.Paired",
    n = nlevels(df[[group_col]])
)
names(fill_cols) <- levels(df[[group_col]])
df$center_light <- .hex_luminance(
    fill_cols[as.character(df[[group_col]])]
) > 0.45

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p <- ggplot(df, aes(x = x, colour = .data[[group_col]])) +
    geom_rect(
        data = bands,
        aes(
            xmin = xmin,
            xmax = xmax,
            ymin = 0,
            ymax = fan_ymax,
            fill = group
        ),
        inherit.aes = FALSE,
        alpha = 0.13,
        colour = NA
    ) +
    geom_rect(
        data = bands,
        aes(
            xmin = xmin,
            xmax = xmax,
            ymin = band_ymin,
            ymax = band_ymax,
            fill = group
        ),
        inherit.aes = FALSE,
        alpha = 0.88,
        colour = NA
    ) +
    geom_segment(
        data = data.frame(x = tick_line_x),
        aes(x = x, xend = x, y = 0, yend = fan_ymax),
        inherit.aes = FALSE,
        colour = "#888888",
        linewidth = 0.55
    ) +
    geom_text(
        data = tick_dat,
        aes(x = x, y = y, label = label),
        inherit.aes = FALSE,
        hjust = -0.15,
        vjust = 0.5,
        colour = "#888888",
        size = 3.8
    ) +
    geom_segment(
        aes(xend = x, y = 0, yend = .data[[mean_col]]),
        linewidth = 0.55
    ) +
    geom_segment(
        aes(
            xend = x,
            y = .data[[mean_col]] - .data[[se_col]],
            yend = .data[[mean_col]] + .data[[se_col]]
        ),
        linewidth = 1.55
    ) +
    geom_segment(
        aes(
            x = x - 0.22,
            xend = x + 0.22,
            y = .data[[mean_col]] - .data[[se_col]],
            yend = .data[[mean_col]] - .data[[se_col]]
        ),
        linewidth = 1.35
    ) +
    geom_segment(
        aes(
            x = x - 0.22,
            xend = x + 0.22,
            y = .data[[mean_col]] + .data[[se_col]],
            yend = .data[[mean_col]] + .data[[se_col]]
        ),
        linewidth = 1.35
    ) +
    geom_point(
        aes(y = .data[[mean_col]], fill = .data[[group_col]]),
        shape = 21,
        size = 5.2,
        stroke = 0.4,
        alpha = 0.96
    ) +
    geom_text(
        data = df[df$center_light, , drop = FALSE],
        aes(
            x = x,
            y = .data[[mean_col]],
            label = .data[[sub_col]]
        ),
        colour = "grey20",
        size = 3.1,
        fontface = "bold",
        hjust = 0.5,
        vjust = 0.5,
        show.legend = FALSE
    ) +
    geom_text(
        data = df[!df$center_light, , drop = FALSE],
        aes(
            x = x,
            y = .data[[mean_col]],
            label = .data[[sub_col]]
        ),
        colour = "white",
        size = 3.1,
        fontface = "bold",
        hjust = 0.5,
        vjust = 0.5,
        show.legend = FALSE
    ) +
    geom_text(
        data = bands,
        aes(
            x = mid,
            y = group_label_y,
            label = group,
            colour = group,
            angle = text_angle
        ),
        inherit.aes = FALSE,
        fontface = "bold",
        size = 6,
        show.legend = FALSE
    ) +
    coord_polar(start = polar_start, direction = 1, clip = "off") +
    scale_x_continuous(limits = c(x_min, x_max), expand = c(0, 0)) +
    scale_y_continuous(
        limits = c(0, group_label_y),
        breaks = y_breaks,
        expand = c(0, 0)
    ) +
    scale_colour_manual(
        values = fill_cols,
        name = config$labels$colour
    ) +
    scale_fill_manual(values = fill_cols, guide = "none") +
    guides(colour = guide_legend(
        override.aes = list(size = 4, linewidth = 1.1),
        order = 2
    )) +
    labs(title = config$labels$title, x = NULL, y = NULL) +
    theme_prism() +
    theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(colour = "#888888", linewidth = 0.4),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        legend.position = "right",
        legend.box = "vertical",
        plot.title = element_text(hjust = 0.5),
        plot.background = element_blank(),
        legend.background = element_blank(),
        plot.margin = margin(24, 65, 24, 65)
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output, width = 9, height = 9)
