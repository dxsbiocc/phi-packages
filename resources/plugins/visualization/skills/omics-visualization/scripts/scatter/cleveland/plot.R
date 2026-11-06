#!/usr/bin/env Rscript

# Template-ID: scatter-cleveland
#
# Purpose:
#   Draw a Cleveland dot plot of a numeric value across categories, with
#   a stem from zero to each point. Optional group colouring, vertical
#   orientation, within-group sorting and values printed in the heads
#   (e.g. GO level-2 term counts by ontology) are CONFIG.
#
# Inputs:
#   A table (tsv/csv). Default example:
#     - Tissue: category label
#     - Count: numeric value
#   Optional group column (CONFIG columns$group), e.g. ontology.
#
# Output:
#   A PDF, PNG, or SVG Cleveland dot plot.
#
# Dependencies:
#   ggplot2, readr
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (columns, group, orientation,
#   sort, head style, value labels, palette, labels). Group colouring
#   and vertical layout are CONFIG, not extra ids.
#
# Scientific assumptions:
#   Each row is one category with a comparable, already-summarised
#   value (not raw replicates). Category order is a display choice.
#   Counts are supplied; this script does not annotate or count genes.

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
        x = "Count",
        y = "Tissue",
        group = NULL             # e.g. "ontology"; NULL = single colour
    ),
    group_order = NULL,          # NULL = first appearance
    orientation = "horizontal",  # "horizontal" or "vertical"
    sort = "input",              # "input" or "desc" (within group)
    head = "ring",               # "ring" or "solid"
    value_labels = FALSE,        # print the value inside each head
    palette = "Qualitative.Bold",
    palette_index = NULL,        # e.g. c(3, 8, 6) to pick hues
    labels = list(
        title = "Cleveland Dot Plot",
        x = "Count",
        y = ""
    ),
    size = list(
        width = 8,
        height = 6
    )
)

load_packages(c("ggplot2", "readr"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- as.data.frame(read_table_auto(io$input), stringsAsFactors = FALSE)
cols <- config$columns
require_columns(df, cols[!vapply(cols, is.null, logical(1))])

x_col <- cols$x
y_col <- cols$y
g_col <- cols$group

df[[x_col]] <- as.numeric(df[[x_col]])
df[[y_col]] <- trimws(as.character(df[[y_col]]))
if (any(!is.finite(df[[x_col]]))) {
    stop(paste("Non-numeric value in", x_col), call. = FALSE)
}

simple <- is.null(g_col) && identical(config$orientation, "horizontal") &&
    identical(config$sort, "input") && identical(config$head, "ring") &&
    !isTRUE(config$value_labels)

if (!simple) {
    if (is.null(g_col)) {
        df$.group <- "all"
        g_col <- ".group"
    } else {
        df[[g_col]] <- trimws(as.character(df[[g_col]]))
    }
    groups <- if (is.null(config$group_order)) unique(df[[g_col]]) else
        as.character(config$group_order)
    df <- df[df[[g_col]] %in% groups, , drop = FALSE]
    groups <- groups[groups %in% df[[g_col]]]
    ord <- if (identical(config$sort, "desc")) {
        order(match(df[[g_col]], groups), -df[[x_col]], seq_len(nrow(df)))
    } else {
        order(match(df[[g_col]], groups), seq_len(nrow(df)))
    }
    df <- df[ord, , drop = FALSE]
    df$pos <- seq_len(nrow(df))
    pal <- palette_colors(config$palette)
    pal <- if (is.null(config$palette_index)) {
        palette_colors(config$palette, n = length(groups))
    } else {
        pal[config$palette_index][seq_along(groups)]
    }
    names(pal) <- groups
    df$col <- unname(pal[df[[g_col]]])
} else {
    df[[y_col]] <- factor(df[[y_col]], levels = unique(df[[y_col]]))
}

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
draw_simple <- function() {
    ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]])) +
        geom_segment(
            aes(x = 0, xend = .data[[x_col]], yend = .data[[y_col]]),
            colour = "#caccd1"
        ) +
        geom_point(size = 5, colour = "#1aafd0", fill = "#ffffffff", shape = 21) +
        geom_point(size = 3, colour = "#1aafd0") +
        labs(
            title = config$labels$title,
            x = config$labels$x,
            y = config$labels$y
        ) +
        theme_minimal(base_size = 14) +
        theme(
            panel.grid.major.y = element_blank(),
            panel.grid.minor = element_blank(),
            axis.title.y = element_blank(),
            axis.title.x = element_text(size = 14, face = "bold"),
            axis.text.y = element_text(size = 12),
            plot.title = element_text(size = 16, face = "bold"),
            plot.background = element_blank()
        )
}

draw_grouped <- function() {
    vertical <- identical(config$orientation, "vertical")
    vmax <- max(df[[x_col]])
    digits <- nchar(format(round(vmax)))
    head_size <- if (isTRUE(config$value_labels)) 4.2 + 1.2 * digits else 4
    text_size <- if (isTRUE(config$value_labels)) 2.9 else 0
    top <- vmax * 1.06
    lo <- -0.04 * top
    lab_at <- lo - 0.012 * top
    breaks <- pretty(c(0, vmax), n = 4)
    breaks <- breaks[breaks <= vmax]

    p <- ggplot(df, aes(x = pos)) +
        geom_segment(aes(xend = pos, y = 0, yend = .data[[x_col]], colour = col),
            linewidth = 1, lineend = "butt", show.legend = FALSE)
    if (identical(config$head, "solid")) {
        p <- p + geom_point(aes(y = .data[[x_col]], fill = .data[[g_col]]),
            shape = 21, colour = "white", stroke = 0, size = head_size)
    } else {
        p <- p +
            geom_point(aes(y = .data[[x_col]], fill = .data[[g_col]]),
                shape = 21, colour = "white", stroke = 0, size = head_size) +
            geom_point(aes(y = .data[[x_col]]), shape = 16, colour = "white",
                size = head_size * 0.55)
    }
    if (isTRUE(config$value_labels)) {
        if (identical(config$head, "solid")) {
            p <- p + geom_text(aes(y = .data[[x_col]], label = .data[[x_col]]),
                colour = "white", fontface = "bold", size = text_size)
        } else {
            p <- p + geom_text(aes(y = .data[[x_col]], label = .data[[x_col]],
                colour = col), fontface = "bold", size = text_size * 0.8)
        }
    }
    p <- p +
        geom_text(aes(y = lab_at, label = .data[[y_col]], colour = col),
            angle = if (vertical) 90 else 0, hjust = 1, vjust = 0.5,
            size = 3.3) +
        scale_colour_identity() +
        scale_fill_manual(values = pal, name = NULL,
            guide = if (identical(g_col, ".group")) "none" else "legend") +
        scale_y_continuous(breaks = breaks) +
        labs(title = config$labels$title, y = config$labels$x, x = NULL) +
        guides(fill = guide_legend(override.aes = list(size = 5))) +
        theme_classic(base_size = 12) +
        theme(
            axis.line = element_blank(),
            axis.ticks = element_line(colour = "black", linewidth = 0.5),
            axis.ticks.length = grid::unit(3, "pt"),
            axis.text.y = element_text(colour = "#4d4d4d", size = 10),
            axis.title.y = element_text(colour = "#4d4d4d", size = 12),
            legend.text = element_text(colour = "#4d4d4d", size = 10),
            legend.key = element_blank(),
            legend.background = element_blank(),
            plot.background = element_blank(),
            plot.title = element_text(size = 14, face = "bold")
        )
    # labels sit outside the panel; reserve room from their real width
    lab_mm <- max(vapply(df[[y_col]], function(s) {
        grid::convertWidth(grid::grobWidth(grid::textGrob(s,
            gp = grid::gpar(fontsize = 3.3 * ggplot2::.pt))), "mm",
            valueOnly = TRUE)
    }, numeric(1)))
    if (vertical) {
        p + scale_x_continuous(breaks = df$pos, labels = NULL,
                expand = expansion(add = 0.7)) +
            coord_cartesian(ylim = c(lo, top), expand = FALSE, clip = "off") +
            theme(plot.margin = margin(10, 10, lab_mm + 6, 10, unit = "mm"),
                axis.ticks.length.x = grid::unit(0, "pt"))
    } else {
        p + coord_flip(ylim = c(lo, top), expand = FALSE, clip = "off") +
            scale_x_reverse(breaks = df$pos, labels = NULL,
                expand = expansion(add = 0.7)) +
            theme(plot.margin = margin(10, 10, 10, lab_mm + 6, unit = "mm"),
                axis.ticks.length.y = grid::unit(0, "pt"))
    }
}

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
p <- if (simple) draw_simple() else draw_grouped()
save_ggplot(p, io$output, width = config$size$width,
    height = config$size$height)
