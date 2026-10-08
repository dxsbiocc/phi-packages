#!/usr/bin/env Rscript

# Template-ID: bar-waffle
#
# Purpose:
#   Draw a faceted waffle chart of integer counts by gene, cancer type,
#   and gender.
#
# Inputs:
#   A table with one row per gene-by-cancer-by-gender. Default example:
#     - Symbol: gene symbol
#     - Cancer_type: cancer type
#     - Gender: sample group
#     - Count: integer count
#
# Output:
#   A PDF, PNG, or SVG waffle chart.
#
# Dependencies:
#   ggplot2, readr
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (column names and labels).
#   Edit DATA PREPARATION to change category order or count rounding.
#   Edit PLOT only when the chart geometry or styling must change.
#
# Scientific assumptions:
#   Each row is a summarized integer count for one gene-by-group cell.
#   Waffle tiles represent counts, not raw replicates.
#   Facets are a display choice for the grouping variables.

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
        fill = "Symbol",
        values = "Count",
        gender = "Gender",
        cancer = "Cancer_type"
    ),
    labels = list(
        title = "Gene mutaion in Pancancer",
        x = "",
        y = ""
    )
)

load_packages(c("ggplot2", "readr"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
require_columns(df, config$columns)

fill_col <- config$columns$fill
value_col <- config$columns$values
gender_col <- config$columns$gender
cancer_col <- config$columns$cancer

df[[fill_col]] <- factor(df[[fill_col]], levels = unique(df[[fill_col]]))
df[[gender_col]] <- factor(df[[gender_col]], levels = unique(df[[gender_col]]))
df[[cancer_col]] <- factor(df[[cancer_col]], levels = unique(df[[cancer_col]]))
df[[value_col]] <- as.integer(round(as.numeric(df[[value_col]])))

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
fill_values <- c(
    "#009392", "#39b185", "#9ccb86", "#e9e29c", "#eeb479", "#e88471", "#cf597e"
)

# geom_waffle(n_rows = 3, flip = TRUE) expands each count into unit cells.
# Within a facet, cells run left to right in rows of 3, then the next row
# up, in data order. The first category sits on the bottom. NA counts are
# dropped, matching na.rm = TRUE.
n_rows <- 3L
panels <- split(
    df,
    interaction(df[[gender_col]], df[[cancer_col]], drop = TRUE)
)
cells <- lapply(panels, function(panel) {
    panel <- panel[!is.na(panel[[value_col]]), , drop = FALSE]
    if (!nrow(panel)) {
        return(NULL)
    }
    index <- rep(seq_len(nrow(panel)), times = panel[[value_col]])
    expanded <- panel[index, , drop = FALSE]
    if (!nrow(expanded)) {
        return(NULL)
    }
    step <- seq_len(nrow(expanded)) - 1L
    expanded$x <- (step %% n_rows) + 1L
    expanded$y <- (step %/% n_rows) + 1L
    expanded
})
cells <- do.call(rbind, cells[!vapply(cells, is.null, logical(1))])

p <- ggplot(cells, aes(
    x = x,
    y = y,
    fill = .data[[fill_col]]
)) +
    geom_tile(color = "white") +
    coord_equal() +
    scale_fill_manual(name = "Genes", values = fill_values) +
    # Cell positions are layout, not data: no breaks, labels, or grid lines.
    scale_x_continuous(breaks = NULL) +
    scale_y_continuous(breaks = NULL) +
    facet_grid(
        rows = vars(.data[[gender_col]]),
        cols = vars(.data[[cancer_col]]),
        switch = "both"
    ) +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y
    ) +
    theme_minimal() +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        strip.background.x = element_blank(),
        # Panels are three cells wide; let the cancer labels overhang them.
        strip.clip = "off",
        strip.background.y = element_rect(fill = "#ffc6c4"),
        panel.background = element_blank(),
        plot.background = element_blank(),
        legend.background = element_blank()
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output)
