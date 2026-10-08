#!/usr/bin/env Rscript

# Template-ID: bar-waterfall
#
# Purpose:
#   Draw a waterfall chart of signed increments from a starting total to an ending total.
#
# Inputs:
#   A two-column table (tsv/csv). Default example:
#     - desc: step label
#     - amount: signed increment; first and last rows are totals
#
# Output:
#   A PDF, PNG, or SVG waterfall chart.
#
# Dependencies:
#   ggplot2, readr, ggprism
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (column names and labels).
#   Edit DATA PREPARATION if first/last rows are not totals.
#   Edit PLOT only when the chart geometry or styling must change.
#
# Scientific assumptions:
#   The first and last rows are totals drawn from zero.
#   Intermediate rows are signed increments applied to a running total.
#   This is a display of already summarized cash-flow steps, not a model.

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
        x = "desc",
        y = "amount"
    ),
    labels = list(
        title = "Waterfall Chart",
        x = NULL,
        y = "Amount"
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
df[[x_col]] <- factor(df[[x_col]], levels = unique(df[[x_col]]))
df[[y_col]] <- as.numeric(df[[y_col]])

n <- nrow(df)
is_total <- seq_len(n) %in% c(1L, n)
running <- 0
ymin <- ymax <- numeric(n)
for (i in seq_len(n)) {
    amt <- df[[y_col]][[i]]
    if (is_total[[i]]) {
        ymin[[i]] <- 0
        ymax[[i]] <- amt
        running <- amt
    } else {
        ymin[[i]] <- running
        ymax[[i]] <- running + amt
        running <- ymax[[i]]
    }
}
df$ymin <- pmin(ymin, ymax)
df$ymax <- pmax(ymin, ymax)
df$kind <- ifelse(is_total, "Total", ifelse(df[[y_col]] >= 0, "Increase", "Decrease"))
df$kind <- factor(df$kind, levels = c("Increase", "Decrease", "Total"))

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p <- ggplot(df, aes(x = .data[[x_col]], fill = kind)) +
    geom_rect(aes(
        xmin = as.numeric(.data[[x_col]]) - 0.4,
        xmax = as.numeric(.data[[x_col]]) + 0.4,
        ymin = ymin,
        ymax = ymax
    )) +
    geom_hline(yintercept = 0, colour = "grey70") +
    geom_text(
        aes(
            y = ymax,
            label = ifelse(
                .data[[y_col]] > 0,
                paste0("+", .data[[y_col]]),
                as.character(.data[[y_col]])
            )
        ),
        vjust = -0.35,
        size = 3
    ) +
    scale_fill_manual(values = c(
        Increase = "#ff6c5f",
        Decrease = "#2dde98",
        Total = "#ff6c5f"
    )) +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y,
        fill = NULL
    ) +
    theme_prism() +
    theme(plot.background = element_blank())

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output)
