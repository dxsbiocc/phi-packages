#!/usr/bin/env Rscript

# Template-ID: bar-venn
#
# Purpose:
#   Draw overlapping Venn circles for 2–4 supplied sets. Counts are
#   unique items in each exclusive region. Five or more sets belong
#   on bar-upset, not a second Venn id.
#
# Inputs:
#   One row per membership. Default example:
#     - item: gene or element id
#     - set: set name
#
# Output:
#   A PDF, PNG, or SVG Venn diagram.
#
# Dependencies:
#   ggplot2, ggvenn, readr, ggprism
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   Extra sets are extra `set` values, not extra ids. 2, 3, or 4
#   circles are CONFIG (the unique names in the table). Percentage
#   labels and listing element names stay in CONFIG. Do not add an
#   id per package (ggvenn / venn). Five or more sets use bar-flower
#   or bar-upset, not a petal Venn on this id.
#
# Scientific assumptions:
#   Membership is supplied. Unique() per set is a display collapse,
#   not a statistical test. Counts are set arithmetic on those ids,
#   not an enrichment p-value.

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
        item = "item",
        set = "set"
    ),
    set_order = c("RNA", "ATAC", "Protein"),
    palettes = list(
        set = "Qualitative.Prism"
    ),
    show_percentage = FALSE,
    show_elements = FALSE,
    labels = list(
        title = ""
    ),
    size = list(
        width = 7.4,
        height = 6.6
    )
)

load_packages(c("ggplot2", "ggvenn", "readr", "ggprism"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- as.data.frame(read_table_auto(io$input), stringsAsFactors = FALSE)
require_columns(df, config$columns)
item_col <- config$columns$item
set_col <- config$columns$set
df[[item_col]] <- as.character(df[[item_col]])
df[[set_col]] <- as.character(df[[set_col]])
ok <- !is.na(df[[item_col]]) & nzchar(trimws(df[[item_col]])) &
    !is.na(df[[set_col]]) & nzchar(trimws(df[[set_col]]))
if (any(!ok)) {
    message("Dropped ", sum(!ok), " row(s) with missing item or set.")
    df <- df[ok, , drop = FALSE]
}
if (!nrow(df)) {
    stop("No set membership left to plot.", call. = FALSE)
}
df[[item_col]] <- trimws(df[[item_col]])
df[[set_col]] <- trimws(df[[set_col]])
df <- unique(df)

sets_in_data <- unique(df[[set_col]])
wanted <- as.character(config$set_order)
wanted <- wanted[wanted %in% sets_in_data]
set_levels <- c(wanted, setdiff(sets_in_data, wanted))
n_set <- length(set_levels)
if (n_set < 2L) {
    stop("A Venn diagram needs at least two sets.", call. = FALSE)
}
if (n_set > 4L) {
    stop(
        paste(
            "This template draws 2–4 overlapping circles.",
            n_set, "sets belong on bar-upset, not a petal Venn."
        ),
        call. = FALSE
    )
}

venn_list <- lapply(set_levels, function(s) {
    unique(df[[item_col]][df[[set_col]] == s])
})
names(venn_list) <- set_levels

pal_all <- palette_colors(config$palettes$set)
# Prism: purple / blue / teal / coral — skip yellow-greens that muddy overlaps.
if (identical(config$palettes$set, "Qualitative.Prism") && length(pal_all) >= 8L) {
    pick <- pal_all[c(1, 2, 3, 8, 9, 10, 11, 4, 5, 6, 7)]
} else {
    pick <- pal_all
}
fill_cols <- expand_palette(pick, n_set)
outline <- palette_colors("Qualitative.Safe")[[12]]

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p <- ggvenn::ggvenn(
    data = venn_list,
    show_elements = isTRUE(config$show_elements),
    show_percentage = isTRUE(config$show_percentage),
    digits = 1,
    fill_color = fill_cols,
    fill_alpha = 0.58,
    stroke_color = outline,
    stroke_size = 0.35,
    set_name_color = "black",
    set_name_size = 6.6,
    text_color = "black",
    text_size = 5
) +
    labs(title = config$labels$title) +
    theme_prism() +
    theme(
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5),
        plot.background = element_blank(),
        panel.background = element_blank(),
        plot.margin = grid::unit(c(6, 8, 6, 8), "mm"),
        legend.position = "none"
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(
    p,
    io$output,
    width = as.numeric(config$size$width),
    height = as.numeric(config$size$height)
)
