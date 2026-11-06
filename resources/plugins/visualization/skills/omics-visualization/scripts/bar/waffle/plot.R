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
#   layout = "parliament" (example_parliament.tsv) is the hemicycle
#   variant: one dot per `unit` counts on concentric arcs, categories
#   filling angular sectors in file order; arc 360 gives a full ring.
#   It needs only columns$fill and columns$values (ggforce).
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
    ),
    # "grid" (default faceted waffle) or "parliament": one dot per unit
    # on concentric arcs, categories filling angular sectors in order.
    # Parliament needs only columns$fill and columns$values.
    layout = "grid",
    parliament = list(
        unit = 1,                    # one dot = this many counts (rounded)
        arc = 180,                   # 180 = hemicycle, 360 = full ring
        rows = NULL,                 # NULL = fewest rows that fit
        inner = 0.4,                 # inner radius / outer radius
        dot_fill = 0.84,             # dot diameter / row spacing
        palette = "Qualitative.Paired / 12",
        border = "black",
        legend_title = NULL,
        size = c(width = 9, height = 5.6)
    )
)

load_packages(c("ggplot2", "readr"))

if (identical(config$layout, "parliament")) {
    pc <- config$parliament
    load_packages("ggforce")
    pdf <- read_table_auto(io$input)
    require_columns(pdf, config$columns[c("fill", "values")])
    cat_lv <- unique(as.character(pdf[[config$columns$fill]]))
    seats <- round(as.numeric(pdf[[config$columns$values]]) / pc$unit)
    if (any(!is.finite(seats) | seats < 0)) {
        stop("values must be finite and non-negative.", call. = FALSE)
    }
    n <- sum(seats)
    if (n < 1) stop("No seats to draw.", call. = FALSE)
    th <- pc$arc / 180 * pi
    fit <- function(r) {
        if (r == 1) return(floor(th / (1 - pc$inner)) + 1)
        d <- (1 - pc$inner) / (r - 1)
        rad <- pc$inner + (seq_len(r) - 1) * d
        sum(floor(th * rad / d) + 1)
    }
    r <- if (is.null(pc$rows)) {
        k <- 1L
        while (fit(k) < n) k <- k + 1L
        k
    } else {
        as.integer(pc$rows)
    }
    rad <- if (r == 1) 1 else seq(pc$inner, 1, length.out = r)
    d <- if (r == 1) 1 - pc$inner else diff(rad)[[1]]
    share <- n * rad / sum(rad)
    per <- floor(share)
    left <- n - sum(per)
    if (left > 0) {
        add <- order(share - per, decreasing = TRUE)[seq_len(left)]
        per[add] <- per[add] + 1
    }
    pos <- do.call(rbind, lapply(seq_len(r), function(i) {
        if (per[[i]] < 1) return(NULL)
        a <- if (per[[i]] == 1) th / 2 else
            seq(0, th, length.out = per[[i]])
        if (pc$arc >= 360) a <- seq(0, th, length.out = per[[i]] + 1)[
            seq_len(per[[i]])]
        data.frame(a = a, rad = rad[[i]])
    }))
    pos <- pos[order(pos$a, -pos$rad), , drop = FALSE]
    pos$cat <- factor(rep(cat_lv, seats), levels = cat_lv)
    pos$x <- pos$rad * cos(pos$a)
    pos$y <- pos$rad * sin(pos$a)
    cols <- stats::setNames(palette_colors(pc$palette, n = length(cat_lv)),
        cat_lv)
    p <- ggplot(pos) +
        ggforce::geom_circle(aes(x0 = x, y0 = y, r = d * pc$dot_fill / 2,
            fill = cat), colour = pc$border, linewidth = 0.4) +
        scale_fill_manual(values = cols, name = if (is.null(pc$legend_title))
            config$columns$fill else pc$legend_title) +
        coord_fixed(clip = "off") +
        labs(title = if (nzchar(config$labels$title)) config$labels$title
            else NULL) +
        theme_void(base_size = 12) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"),
            plot.background = element_blank(),
            legend.key.size = unit(6, "mm"))
    save_ggplot(p, io$output, width = pc$size[["width"]],
        height = pc$size[["height"]])
    quit(save = "no", status = 0)
}

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
