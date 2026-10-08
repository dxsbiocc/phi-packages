#!/usr/bin/env Rscript

# Template-ID: scatter-dumbbell-region
#
# Purpose:
#   Draw a dumbbell plot of two groups on a shared categorical axis, with
#   shaded mean +/- SD regions behind the points.
#
# Inputs:
#   A table with one row per category-by-group. Default example:
#     - age_group: category label
#     - mean_prev: numeric value
#     - Sex: group (two levels)
#
# Output:
#   A PDF, PNG, or SVG regional dumbbell plot.
#
# Dependencies:
#   ggplot2, readr, ggprism
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (column names and labels).
#   Edit DATA PREPARATION to change how mean and SD regions are
#   computed, or category order.
#   Edit PLOT only when the chart geometry or styling must change.
#
# Scientific assumptions:
#   Each row is one category-by-group pair with a comparable numeric value.
#   Mean and SD bands summarize the plotted values, not a formal test.
#   Category order is a display choice.

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
        x = "mean_prev",
        y = "age_group",
        group = "Sex"
    ),
    labels = list(
        title = "Mean Prevalence of obesity",
        x = "",
        y = ""
    )
)

load_packages(c("ggplot2", "readr", "ggprism"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
require_columns(df, config$columns)

x_col <- config$columns$x
y_col <- config$columns$y
group_col <- config$columns$group

df[[y_col]] <- factor(df[[y_col]], levels = unique(df[[y_col]]))
df[[group_col]] <- factor(df[[group_col]], levels = unique(df[[group_col]]))
df[[x_col]] <- as.numeric(df[[x_col]])

summary_data <- do.call(rbind, lapply(
    split(df, df[[group_col]]),
    function(sub) {
        mu <- mean(sub[[x_col]], na.rm = TRUE)
        se <- stats::sd(sub[[x_col]], na.rm = TRUE)
        data.frame(
            mu = mu,
            se = se,
            upper = mu + se,
            lower = mu - se,
            group = sub[[group_col]][[1]]
        )
    }
))
names(summary_data)[names(summary_data) == "group"] <- group_col
summary_data[[group_col]] <- factor(
    summary_data[[group_col]],
    levels = levels(df[[group_col]])
)

last_y <- as.character(tail(levels(df[[y_col]]), 1))
region_idx <- nrow(summary_data)

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
colour_values <- c("#f784b6", "#79ceb8")
fill_values <- c("#FFDEEA", "#E5FAED")

p <- ggplot(df, aes(
    x = .data[[config$columns$x]],
    y = .data[[config$columns$y]],
    colour = .data[[config$columns$group]]
)) +
    geom_rect(
        data = summary_data,
        aes(
            xmin = lower,
            xmax = upper,
            fill = .data[[config$columns$group]]
        ),
        ymin = -Inf,
        ymax = Inf,
        alpha = 0.2,
        inherit.aes = FALSE
    ) +
    annotate(
        "text",
        x = summary_data$mu[[region_idx]],
        y = last_y,
        label = "MEAN",
        angle = 90,
        size = 3,
        vjust = -0.2,
        colour = "#79ceb8"
    ) +
    annotate(
        "text",
        x = summary_data$upper[[region_idx]],
        y = last_y,
        label = "STDEV",
        angle = 90,
        size = 3,
        vjust = -0.2,
        colour = "#79ceb8"
    ) +
    geom_vline(
        data = summary_data,
        aes(
            xintercept = mu,
            colour = .data[[config$columns$group]]
        ),
        show.legend = FALSE
    ) +
    geom_line(
        aes(group = .data[[config$columns$y]]),
        colour = "#cba7a2"
    ) +
    geom_point(size = 5, shape = 21, fill = "#ffffff") +
    geom_point(size = 3) +
    scale_colour_manual(values = colour_values) +
    scale_fill_manual(values = fill_values) +
    scale_y_discrete(expand = expansion(mult = c(0.1, 0.1))) +
    guides(fill = "none") +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y
    ) +
    theme_prism() +
    theme(
        axis.line = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        plot.background = element_blank(),
        legend.position = "top"
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output)
