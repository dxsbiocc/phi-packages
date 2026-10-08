#!/usr/bin/env Rscript

# Template-ID: scatter-jitter-group
#
# Purpose:
#   Draw a grouped jitter scatter plot of log2 fold-change across
#   contrasts, with category tiles and labels for extreme genes.
#
# Inputs:
#   A table with one row per gene-by-contrast. Default example:
#     - cate: contrast label
#     - logFC: log2 fold-change
#     - sign: direction (Up / Down)
#     - symbol: gene label
#
# Output:
#   A PDF, PNG, or SVG grouped jitter plot.
#
# Dependencies:
#   ggplot2, readr, ggprism, ggrepel, ggnewscale
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (column names and labels).
#   Edit DATA PREPARATION to change contrast order, jitter seed, or
#   how many labels are drawn per group.
#   Edit PLOT only when the chart geometry or styling must change.
#
# Scientific assumptions:
#   Each row is one gene in one contrast with a log2 fold-change.
#   Jitter is a display transform to reduce overplotting.
#   Background bars use the min/max of each direction, not a test.

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
        x = "cate",
        y = "logFC",
        group = "sign",
        label = "symbol"
    ),
    labels = list(
        title = "Yap1 Overexpression and Knockout",
        x = "",
        y = "log2FC"
    )
)

load_packages(c("ggplot2", "readr", "ggprism", "ggrepel", "ggnewscale"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
require_columns(df, config$columns)

x_col <- config$columns$x
y_col <- config$columns$y
group_col <- config$columns$group
label_col <- config$columns$label

keep <- c(x_col, y_col, group_col, label_col)
df <- df[stats::complete.cases(df[, keep]), ]
df[[x_col]] <- factor(df[[x_col]], levels = unique(df[[x_col]]))
df[[group_col]] <- factor(df[[group_col]], levels = unique(df[[group_col]]))
df[[y_col]] <- as.numeric(df[[y_col]])

df$hight <- NA_real_
for (g in levels(df[[group_col]])) {
    idx <- df[[group_col]] == g
    v <- df[[y_col]][idx]
    df$hight[idx] <- ifelse(df[[y_col]][idx] < 0, min(v), max(v))
}

col_data <- unique(df[, c(x_col, "hight"), drop = FALSE])
tile_data <- unique(df[, x_col, drop = FALSE])

p_tmp <- ggplot(df, aes(
    x = .data[[x_col]],
    y = .data[[y_col]]
)) +
    geom_point(position = position_jitter(width = 0.4, height = 0, seed = 1234))
built <- ggplot_build(p_tmp)$data[[1]]
point_data <- df
point_data[[x_col]] <- built$x
point_data[[y_col]] <- built$y

text_data <- do.call(rbind, lapply(
    split(
        point_data,
        list(point_data[[group_col]], point_data[[y_col]] > 0),
        drop = TRUE
    ),
    function(sub) {
        sub[head(order(abs(sub[[y_col]]), decreasing = TRUE), 5), , drop = FALSE]
    }
))

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
sign_fill <- c(Down = "#2dde98", Up = "#ff6c5f")
cate_fill <- c("#3be8b0", "#1aafd0", "#6a67ce", "#ffb900", "#fc636b")
cate_colour <- c("#595d7a", "#595d7a", "#c5c5c5", "#595d7a", "#595d7aff")

p <- ggplot(df, aes(
    x = .data[[config$columns$x]],
    y = .data[[config$columns$y]],
    colour = .data[[config$columns$group]]
)) +
    geom_col(
        data = col_data,
        aes(x = .data[[config$columns$x]], y = hight),
        fill = "#ebf3f5ff",
        alpha = 0.6,
        width = 0.8,
        inherit.aes = FALSE,
        show.legend = FALSE
    ) +
    geom_point(
        data = point_data,
        aes(fill = .data[[config$columns$group]], group = .data[[config$columns$x]]),
        colour = "#8a8b8c",
        shape = 21,
        stroke = 0
    ) +
    scale_fill_manual(
        name = "Change",
        values = sign_fill,
        labels = c("Down", "Up"),
        drop = FALSE,
        guide = guide_legend(override.aes = list(size = 4))
    ) +
    new_scale_fill() +
    geom_tile(
        data = tile_data,
        aes(x = .data[[config$columns$x]], y = 0, fill = .data[[config$columns$x]]),
        height = 1.4,
        colour = "black",
        alpha = 0.6,
        width = 0.8,
        show.legend = FALSE
    ) +
    geom_text(
        data = tile_data,
        aes(
            x = .data[[config$columns$x]],
            y = 0,
            label = .data[[config$columns$x]],
            colour = .data[[config$columns$x]]
        ),
        size = 3.5,
        show.legend = FALSE,
        inherit.aes = FALSE
    ) +
    geom_text_repel(
        data = text_data,
        aes(
            x = .data[[config$columns$x]],
            y = .data[[config$columns$y]],
            label = .data[[config$columns$label]]
        ),
        force = 2,
        size = 3,
        arrow = arrow(length = unit(0.008, "npc"), type = "open", ends = "last"),
        inherit.aes = FALSE
    ) +
    scale_fill_manual(values = cate_fill) +
    scale_colour_manual(values = cate_colour) +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y
    ) +
    theme_prism() +
    theme(
        plot.background = element_blank(),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 13, colour = "black", face = "bold"),
        axis.line.y = element_line(colour = "black", linewidth = 1),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank()
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output)
