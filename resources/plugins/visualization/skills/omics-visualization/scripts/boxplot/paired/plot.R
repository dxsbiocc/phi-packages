#!/usr/bin/env Rscript

# Template-ID: boxplot-paired
#
# Purpose:
#   Draw paired boxplots with Bezier curves connecting matched
#   observations and a paired Wilcoxon significance bracket.
#
# Inputs:
#   A table with one row per sample. Default example:
#     - patient: pairing id
#     - sample_type: category (two levels)
#     - BRCA1: numeric value
#
# Output:
#   A PDF, PNG, or SVG paired boxplot.
#
# Dependencies:
#   ggplot2, readr, ggforce, ggprism
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (column names and labels).
#   Edit DATA PREPARATION to change pairing or category order.
#   Edit PLOT only when the chart geometry or styling must change.
#
# Scientific assumptions:
#   Rows that share the same pairing id are matched observations.
#   Bezier curves are a display transform connecting paired values.
#   The bracket is a paired Wilcoxon rank-sum test at p < 0.05.

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
        x = "sample_type",
        y = "BRCA1",
        id = "patient"
    ),
    labels = list(
        title = "TCGA-ACC",
        x = "",
        y = "Expression of BRCA1"
    )
)

load_packages(c("ggplot2", "readr", "ggforce", "ggprism"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
require_columns(df, config$columns)
x_col <- config$columns$x
y_col <- config$columns$y
id_col <- config$columns$id

df[[x_col]] <- factor(df[[x_col]], levels = unique(df[[x_col]]))
df[[y_col]] <- as.numeric(df[[y_col]])
df$x_plot <- as.numeric(df[[x_col]])
x_lvls <- levels(df[[x_col]])

# Bezier control points connecting paired observations (by id).
bezier_rows <- list()
k <- 1L
for (pid in unique(df[[id_col]])) {
    sub <- df[df[[id_col]] == pid, ]
    sub <- sub[order(sub$x_plot), ]
    if (nrow(sub) < 2L) {
        next
    }
    for (i in seq_len(nrow(sub) - 1L)) {
        bezier_rows[[k]] <- data.frame(
            bx = c(
                sub$x_plot[[i]], sub$x_plot[[i]] + 0.3,
                sub$x_plot[[i + 1L]] - 0.3, sub$x_plot[[i + 1L]]
            ),
            by = c(
                sub[[y_col]][[i]], sub[[y_col]][[i]],
                sub[[y_col]][[i + 1L]], sub[[y_col]][[i + 1L]]
            ),
            bgroup = paste0(pid, "_", i)
        )
        k <- k + 1L
    }
}
bezier_data <- do.call(rbind, bezier_rows)

comparisons <- list()
for (i in seq_len(length(x_lvls) - 1L)) {
    for (j in seq(i + 1L, length(x_lvls))) {
        left <- df[df[[x_col]] == x_lvls[[i]], c(id_col, y_col)]
        right <- df[df[[x_col]] == x_lvls[[j]], c(id_col, y_col)]
        names(left) <- c("id", "a")
        names(right) <- c("id", "b")
        merged <- merge(left, right, by = "id")
        if (nrow(merged) < 1L) {
            next
        }
        pv <- tryCatch(
            wilcox.test(merged$a, merged$b, paired = TRUE)$p.value,
            error = function(e) 1
        )
        if (is.finite(pv) && pv < 0.05) {
            comparisons[[length(comparisons) + 1L]] <- list(
                xmin = i, xmax = j, pvalue = pv
            )
        }
    }
}

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
fill_values <- c("#69c5e4", "#e66760")

p <- ggplot(df, aes(
    x = x_plot,
    y = .data[[y_col]]
)) +
    geom_bezier(
        data = bezier_data,
        aes(x = bx, y = by, group = bgroup),
        colour = "#caccd1",
        alpha = 0.5,
        inherit.aes = FALSE
    ) +
    geom_boxplot(
        aes(colour = .data[[x_col]]),
        fill = "#ffffff",
        alpha = 0.11,
        outlier.shape = NA,
        show.legend = FALSE
    ) +
    geom_point(
        position = position_dodge2(width = 0.3),
        size = 4,
        colour = "white"
    ) +
    geom_point(
        aes(fill = .data[[x_col]]),
        position = position_dodge2(width = 0.3),
        size = 3,
        shape = 21,
        stroke = 0.5,
        alpha = 0.9,
        show.legend = FALSE
    ) +
    scale_fill_manual(values = fill_values) +
    scale_colour_manual(values = fill_values) +
    scale_x_continuous(
        breaks = seq_along(x_lvls),
        labels = x_lvls
    ) +
    scale_y_continuous(expand = c(0.1, 0.1)) +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y
    ) +
    theme_prism() +
    theme(
        plot.background = element_blank(),
        plot.title = element_text(size = 14, hjust = 0.5)
    )

if (length(comparisons) > 0) {
    y_max <- max(df[[y_col]], na.rm = TRUE)
    y_span <- diff(range(df[[y_col]], na.rm = TRUE))
    for (i in seq_along(comparisons)) {
        cmp <- comparisons[[i]]
        y_pos <- y_max + 0.12 * y_span * i
        tip <- 0.03 * y_span
        p <- p +
            annotate(
                "segment",
                x = cmp$xmin, xend = cmp$xmax,
                y = y_pos, yend = y_pos
            ) +
            annotate(
                "segment",
                x = cmp$xmin, xend = cmp$xmin,
                y = y_pos - tip, yend = y_pos
            ) +
            annotate(
                "segment",
                x = cmp$xmax, xend = cmp$xmax,
                y = y_pos - tip, yend = y_pos
            ) +
            annotate(
                "text",
                x = mean(c(cmp$xmin, cmp$xmax)),
                y = y_pos,
                label = sprintf("pvalue = %.2g", cmp$pvalue),
                vjust = -0.4,
                size = 3.5
            )
    }
}

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output)
