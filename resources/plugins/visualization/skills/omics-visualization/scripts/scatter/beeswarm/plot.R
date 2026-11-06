#!/usr/bin/env Rscript

# Template-ID: scatter-beeswarm
#
# Purpose:
#   Draw a beeswarm / quasirandom plot of a numeric value across
#   categories. style = "quartile" (default) overlays IQR bars and a
#   median marker, optional second colour group. style = "mean_sem"
#   overlays mean ± SEM as a short bar + whiskers, optional supplied
#   significance brackets, and a white classic panel (preview_ccl2.png).
#
# Inputs:
#   One row per observation. Default example:
#     - name, value, group
#   mean_sem example (example_ccl2.tsv):
#     - group, value
#   Optional sidecar for mean_sem:
#     - brackets.tsv: group1, group2, label
#
# Output:
#   A PDF, PNG, or SVG beeswarm plot.
#
# Dependencies:
#   ggplot2, readr, ggprism, ggbeeswarm, scales
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   Edit CONFIG for columns, colours, summary style, brackets, and
#   labels. Extra categories are extra levels, not a new id. A box
#   plus swarm stays boxplot-basic. Random jitter without beeswarm
#   layout stays scatter-jitter. Do not run ANOVA or Tukey in this
#   script; bracket labels are supplied.
#
# Scientific assumptions:
#   Each row is one observation. Quartiles / mean / SEM are display
#   summaries of the supplied values, not a hypothesis test.

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
    # "quartile": IQR bar + median marker (default preview).
    # "mean_sem": mean ± SEM bar + optional brackets (preview_ccl2.png).
    style = "quartile",
    columns = list(
        x = "name",
        y = "value",
        # Optional second colour group for quartile style; NULL = colour by x.
        group = "group"
    ),
    group_order = NULL,
    palettes = list(group = "Qualitative.Safe"),
    # Named colours keyed by x (or group) level; overrides palette.
    group_colours = NULL,
    point_size = 2.8,
    point_alpha = 0.75,
    point_stroke = 0.25,
    swarm_width = 0.18,
    mean_sem = list(
        mean_width = 0.32,
        mean_linewidth = 1.1,
        error_linewidth = 0.55,
        error_width = 0.1,
        brackets_file = "brackets_ccl2.tsv",
        bracket_pad = 2.2,
        bracket_tick = 0.7,
        bracket_label_size = 5.5,
        panel_label = "i",
        x_angle = 45,
        width = 4.0,
        height = 4.6
    ),
    labels = list(
        title = NULL,
        x = NULL,
        y = "Expression level (log2 TPM)",
        colour = NULL
    ),
    size = list(width = 7, height = 5)
)

load_packages(c("ggplot2", "readr", "ggprism", "ggbeeswarm", "scales"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- as.data.frame(read_table_auto(io$input), stringsAsFactors = FALSE)
x_col <- config$columns$x
y_col <- config$columns$y
g_col <- config$columns$group
require_columns(df, list(x = x_col, y = y_col))

df[[x_col]] <- as.character(df[[x_col]])
df[[y_col]] <- as.numeric(df[[y_col]])
df <- df[is.finite(df[[y_col]]) & !is.na(df[[x_col]]) & nzchar(df[[x_col]]), , drop = FALSE]
if (!nrow(df)) stop("No finite rows to plot.", call. = FALSE)

x_lv <- if (!is.null(config$group_order)) {
    c(
        intersect(as.character(config$group_order), unique(df[[x_col]])),
        setdiff(unique(df[[x_col]]), as.character(config$group_order))
    )
} else {
    unique(df[[x_col]])
}
df[[x_col]] <- factor(df[[x_col]], levels = x_lv)

use_group <- !is.null(g_col) && g_col %in% names(df) &&
    !identical(g_col, x_col)
if (use_group) {
    df[[g_col]] <- factor(as.character(df[[g_col]]), levels = unique(as.character(df[[g_col]])))
    colour_col <- g_col
} else {
    colour_col <- x_col
}

col_lv <- levels(df[[colour_col]])
if (!is.null(config$group_colours) && length(config$group_colours)) {
    g_cols <- unlist(config$group_colours)
    miss <- setdiff(col_lv, names(g_cols))
    if (length(miss)) {
        extra <- expand_palette(palette_colors(config$palettes$group), length(miss))
        names(extra) <- miss
        g_cols <- c(g_cols, extra)
    }
    g_cols <- g_cols[col_lv]
} else {
    g_cols <- expand_palette(palette_colors(config$palettes$group), length(col_lv))
    names(g_cols) <- col_lv
}

sidecar_path <- function(name) {
    if (is.null(name) || !nzchar(as.character(name))) return(NULL)
    path <- if (file.exists(name)) name else file.path(dirname(io$input), name)
    if (file.exists(path)) path else NULL
}

# -----------------------------------------------------------------------------
# PLOT — mean ± SEM
# -----------------------------------------------------------------------------
if (identical(config$style, "mean_sem")) {
    ms <- config$mean_sem
    sum_df <- do.call(rbind, lapply(split(df, df[[x_col]], drop = TRUE), function(sub) {
        v <- sub[[y_col]]
        n <- length(v)
        m <- mean(v)
        se <- if (n > 1L) stats::sd(v) / sqrt(n) else 0
        data.frame(
            x = sub[[x_col]][[1]],
            mean = m,
            ymin = m - se,
            ymax = m + se,
            stringsAsFactors = FALSE
        )
    }))
    sum_df$x <- factor(as.character(sum_df$x), levels = x_lv)

    p <- ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]])) +
        ggbeeswarm::geom_quasirandom(
            aes(fill = .data[[colour_col]]),
            width = config$swarm_width,
            size = config$point_size,
            shape = 21,
            colour = "grey20",
            stroke = config$point_stroke,
            alpha = config$point_alpha,
            show.legend = FALSE
        ) +
        geom_errorbar(
            data = sum_df,
            aes(x = x, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE,
            width = ms$error_width,
            linewidth = ms$error_linewidth,
            colour = "grey15"
        ) +
        geom_segment(
            data = sum_df,
            aes(
                x = as.numeric(x) - ms$mean_width / 2,
                xend = as.numeric(x) + ms$mean_width / 2,
                y = mean, yend = mean
            ),
            inherit.aes = FALSE,
            linewidth = ms$mean_linewidth,
            colour = "grey10"
        ) +
        scale_fill_manual(values = g_cols, guide = "none") +
        scale_x_discrete(drop = FALSE) +
        scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
        labs(
            title = config$labels$title,
            x = config$labels$x,
            y = config$labels$y
        ) +
        theme_classic(base_size = 12) +
        theme(
            axis.text = element_text(colour = "black"),
            axis.title = element_text(colour = "black"),
            axis.text.x = element_text(
                angle = ms$x_angle,
                hjust = if (ms$x_angle > 0) 1 else 0.5,
                vjust = 1
            ),
            axis.ticks = element_line(colour = "black"),
            panel.grid = element_blank(),
            panel.background = element_rect(fill = "white", colour = NA),
            plot.background = element_rect(fill = "white", colour = NA),
            plot.margin = margin(10, 10, 8, 10)
        )

    if (!is.null(ms$panel_label) && nzchar(as.character(ms$panel_label))) {
        p <- p +
            annotate(
                "text",
                x = -Inf, y = Inf,
                label = ms$panel_label,
                hjust = -0.4, vjust = 1.4,
                size = 5, fontface = "bold", colour = "grey10"
            )
    }

    br_path <- sidecar_path(ms$brackets_file)
    if (!is.null(br_path)) {
        br <- as.data.frame(read_table_auto(br_path), stringsAsFactors = FALSE)
        require_columns(br, list("group1", "group2", "label"))
        br$x1 <- match(as.character(br$group1), x_lv)
        br$x2 <- match(as.character(br$group2), x_lv)
        if (anyNA(br$x1) || anyNA(br$x2)) {
            stop("brackets.tsv groups must match x levels.", call. = FALSE)
        }
        y_top <- max(c(df[[y_col]], sum_df$ymax), na.rm = TRUE)
        br$y <- y_top + ms$bracket_pad +
            (seq_len(nrow(br)) - 1L) * (ms$bracket_tick * 3)
        tk <- ms$bracket_tick
        seg <- rbind(
            data.frame(x = br$x1, xend = br$x2, y = br$y, yend = br$y),
            data.frame(x = br$x1, xend = br$x1, y = br$y - tk, yend = br$y),
            data.frame(x = br$x2, xend = br$x2, y = br$y - tk, yend = br$y)
        )
        p <- p +
            geom_segment(
                data = seg,
                aes(x = x, xend = xend, y = y, yend = yend),
                inherit.aes = FALSE,
                linewidth = 0.4, colour = "grey20"
            ) +
            geom_text(
                data = br,
                aes(x = (x1 + x2) / 2, y = y + tk * 0.35, label = label),
                inherit.aes = FALSE,
                size = ms$bracket_label_size,
                colour = "grey10",
                vjust = 0,
                fontface = "bold"
            ) +
            coord_cartesian(
                ylim = c(0, max(br$y) + tk * 2.2),
                clip = "off"
            )
    } else {
        p <- p + coord_cartesian(ylim = c(0, NA), clip = "off")
    }

    if (is.null(config$labels$x) || !nzchar(as.character(config$labels$x))) {
        p <- p + theme(axis.title.x = element_blank())
    }
    if (is.null(config$labels$title) || !nzchar(as.character(config$labels$title))) {
        p <- p + theme(plot.title = element_blank())
    }

    save_ggplot(p, io$output, width = ms$width, height = ms$height)
    quit(save = "no")
}

# -----------------------------------------------------------------------------
# PLOT — quartile (legacy default)
# -----------------------------------------------------------------------------
if (!use_group) {
    # Mirror legacy API: require a group column when using quartile dodge.
    if (is.null(g_col) || !g_col %in% names(df)) {
        df$.group <- df[[x_col]]
        g_col <- ".group"
        colour_col <- g_col
        df[[g_col]] <- df[[x_col]]
        g_cols <- expand_palette(palette_colors(config$palettes$group), nlevels(df[[x_col]]))
        names(g_cols) <- levels(df[[x_col]])
    }
}

summary_data <- do.call(rbind, lapply(
    split(df, list(df[[x_col]], df[[g_col]]), drop = TRUE),
    function(sub) {
        v <- sub[[y_col]]
        data.frame(
            name = sub[[x_col]][[1]],
            group = sub[[g_col]][[1]],
            first = as.numeric(stats::quantile(v, 0.25)),
            second = as.numeric(stats::quantile(v, 0.5)),
            third = as.numeric(stats::quantile(v, 0.75))
        )
    }
))
names(summary_data)[1:2] <- c(x_col, g_col)
summary_data[[x_col]] <- factor(summary_data[[x_col]], levels = levels(df[[x_col]]))
summary_data[[g_col]] <- factor(summary_data[[g_col]], levels = levels(df[[g_col]]))

p <- ggplot(df, aes(
    x = .data[[x_col]],
    y = .data[[y_col]],
    colour = .data[[g_col]]
)) +
    geom_quasirandom(width = 0.15, size = 3, dodge.width = 0.8, alpha = 0.6) +
    geom_errorbar(
        data = summary_data,
        aes(
            x = .data[[x_col]],
            ymin = first,
            ymax = third,
            group = .data[[g_col]]
        ),
        width = 0,
        position = position_dodge(width = 0.8),
        linewidth = 1,
        inherit.aes = FALSE
    ) +
    geom_point(
        data = summary_data,
        aes(
            x = .data[[x_col]],
            y = second,
            group = .data[[g_col]]
        ),
        position = position_dodge2(width = 0.8),
        size = 5,
        shape = 21,
        fill = "white",
        stroke = 1,
        inherit.aes = FALSE
    ) +
    geom_point(
        data = summary_data,
        aes(
            x = .data[[x_col]],
            y = second,
            fill = .data[[g_col]],
            group = .data[[g_col]]
        ),
        position = position_dodge2(width = 0.8),
        size = 4,
        shape = 21,
        stroke = 0,
        inherit.aes = FALSE
    ) +
    scale_colour_manual(values = g_cols) +
    scale_fill_manual(values = g_cols) +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y,
        colour = config$labels$colour,
        fill = config$labels$colour
    ) +
    theme_prism() +
    theme(
        legend.position = "inside",
        legend.position.inside = c(0.9, 0.8),
        legend.background = element_rect(fill = "white", colour = "grey80"),
        axis.line = element_blank(),
        panel.border = element_rect(fill = "transparent", colour = "black"),
        plot.background = element_blank()
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(
    p,
    io$output,
    width = config$size$width,
    height = config$size$height
)
