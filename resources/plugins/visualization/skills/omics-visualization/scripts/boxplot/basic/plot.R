#!/usr/bin/env Rscript

# Template-ID: boxplot-basic
#
# Purpose:
#   Draw a basic boxplot comparing a numeric value across categories,
#   with jittered points overlaid. style = "beeswarm" draws outline
#   boxes, a quasirandom swarm, optional log10 y, multi-line x labels
#   from a sidecar, and supplied comparison brackets. style = "notch"
#   draws notched boxes (median CI display), coloured outlines / points,
#   optional horizontal orientation (preview_notch.png). style =
#   "cohort" draws treatment × day boxes with day-coloured fill/points,
#   dotted cohort dividers, underlined treatment labels, optional log10
#   y, and supplied brackets (preview_cutrun.png).
#
# Inputs:
#   A two-column table (tsv/csv). Default example:
#     - Country: category label
#     - Income: numeric value
#   beeswarm example (example_beeswarm.tsv):
#     - group, value
#   notch example (example_notch.tsv):
#     - cell_type, ratio
#   cohort example (example_cutrun.tsv):
#     - treatment, day, reads
#   Optional sidecars next to the input:
#     - labels.tsv: group, label  (| becomes a line break)
#     - brackets.tsv: group1, group2, label  (fold / P text is supplied)
#       For cohort style, group keys are "treatment|day".
#
# Output:
#   A PDF, PNG, or SVG boxplot.
#
# Dependencies:
#   ggplot2, readr, ggprism, scales; beeswarm also needs ggbeeswarm
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (column names and labels).
#   Edit DATA PREPARATION to change category order.
#   Edit PLOT only when the chart geometry or styling must change.
#   Do not run fold-change or significance tests in this script;
#   bracket labels come from brackets.tsv. notch / cohort layout are
#   CONFIG; do not add an id per assay or treatment panel.
#
# Scientific assumptions:
#   Each row is one observation of a numeric value in a category.
#   Boxes summarize the distribution; points are the raw values.
#   Notches are a display of the median confidence interval from
#   geom_boxplot, not a hypothesis test run in this script.

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
    # "basic": filled boxes + jitter (default preview).
    # "beeswarm": outline boxes + quasirandom points, white panel,
    #            optional log10 y and supplied brackets.
    # "notch": notched boxes, coloured outline + jitter, white panel;
    #          orientation "h" matches immune-cell frequency layouts.
    # "cohort": treatment × day boxes, day colours, cohort dividers,
    #           underlined treatment labels, optional brackets.
    style = "basic",
    columns = list(
        x = "Country",
        y = "Income",
        # cohort style (ignored otherwise):
        treatment = "treatment",
        day = "day"
    ),
    labels = list(
        title = "Basic_Boxplot",
        x = "Country",
        y = "Income"
    ),
    group_order = NULL,
    cohort = list(
        treatment_order = c("PBT", "PBT+IFNa2a", "IFNa2a"),
        day_order = c("D0", "D4", "D28"),
        day_colours = c(
            D0 = "#7CB342",
            D4 = "#8E3A59",
            D28 = "#E07A3D"
        ),
        palette = "Qualitative.Bold",
        y_trans = "log10",
        box_alpha = 0.35,
        box_width = 0.55,
        point_size = 1.35,
        point_alpha = 0.75,
        jitter_width = 0.18,
        divider_linetype = "dotted",
        divider_colour = "grey40",
        brackets_file = "brackets_cutrun.tsv",
        bracket_tick = 0.04,
        bracket_pad = 0.12,
        bracket_step = 0.22,
        label_size = 4.2,
        treatment_label_size = 3.4,
        base_size = 12,
        width = 8.5,
        height = 6.2
    ),
    beeswarm = list(
        palette = "Qualitative.Pastel",
        # One catalog colour for the median line.
        median = "Qualitative.Bold:5",
        box_colour = "grey25",
        box_linewidth = 0.55,
        point_size = 4.2,
        point_alpha = 0.9,
        point_width = 0.16,
        y_trans = "log10",  # "log10" | "identity"
        labels_file = "labels.tsv",
        brackets_file = "brackets.tsv",
        bracket_tick = 0.04,
        bracket_pad = 0.08,  # log10 (or raw) pad above the tallest point
        bracket_step = 0.22, # vertical gap between stacked brackets
        label_size = 3.2,
        base_size = 13,
        grid = TRUE,
        width = 8.5,
        height = 6.2
    ),
    notch = list(
        palette = "Qualitative.Bold",
        # "h": category on y, value on x. "v": category on x, value on y.
        orientation = "h",
        notch = TRUE,
        notchwidth = 0.5,
        fill = "white",          # box fill; outline uses category colour
        fill_alpha = 1,
        box_linewidth = 0.55,
        outlier.shape = NA,      # raw points shown via jitter
        point_size = 1.6,
        point_alpha = 0.65,
        jitter_width = 0.15,     # along the categorical axis
        base_size = 12,
        width = 8.5,
        height = 9.0
    )
)

load_packages(c("ggplot2", "readr", "ggprism", "scales"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- as.data.frame(read_table_auto(io$input), stringsAsFactors = FALSE)

pal_ref <- function(ref) {
    parts <- strsplit(ref, ":", fixed = TRUE)[[1]]
    cols <- palette_colors(parts[[1]])
    if (length(parts) > 1L) cols[[as.integer(parts[[2]])]] else cols
}

sidecar_path <- function(name) {
    if (is.null(name) || !nzchar(as.character(name))) return(NULL)
    path <- if (file.exists(name)) name else file.path(dirname(io$input), name)
    if (file.exists(path)) path else NULL
}

# -----------------------------------------------------------------------------
# PLOT — cohort (treatment × day)
# -----------------------------------------------------------------------------
if (identical(config$style, "cohort")) {
    tr_col <- config$columns$treatment
    day_col <- config$columns$day
    y_col <- config$columns$y
    require_columns(df, list(treatment = tr_col, day = day_col, y = y_col))
    cc <- config$cohort
    df[[y_col]] <- as.numeric(df[[y_col]])
    df[[tr_col]] <- as.character(df[[tr_col]])
    df[[day_col]] <- as.character(df[[day_col]])
    df <- df[is.finite(df[[y_col]]), , drop = FALSE]
    if (!nrow(df)) stop("No finite rows to plot.", call. = FALSE)

    tr_lv <- if (!is.null(cc$treatment_order)) {
        c(
            intersect(as.character(cc$treatment_order), unique(df[[tr_col]])),
            setdiff(unique(df[[tr_col]]), as.character(cc$treatment_order))
        )
    } else {
        unique(df[[tr_col]])
    }
    day_lv <- if (!is.null(cc$day_order)) {
        c(
            intersect(as.character(cc$day_order), unique(df[[day_col]])),
            setdiff(unique(df[[day_col]]), as.character(cc$day_order))
        )
    } else {
        unique(df[[day_col]])
    }
    df[[tr_col]] <- factor(df[[tr_col]], levels = tr_lv)
    df[[day_col]] <- factor(df[[day_col]], levels = day_lv)
    # Combined x key: treatment|day in treatment-major order.
    x_lv <- as.vector(t(outer(tr_lv, day_lv, paste, sep = "|")))
    df$.x <- factor(
        paste(as.character(df[[tr_col]]), as.character(df[[day_col]]), sep = "|"),
        levels = x_lv
    )
    df <- df[!is.na(df$.x), , drop = FALSE]

    day_cols <- cc$day_colours
    if (is.null(day_cols)) day_cols <- character()
    day_cols <- unlist(day_cols)
    miss_d <- setdiff(day_lv, names(day_cols))
    if (length(miss_d)) {
        extra <- expand_palette(palette_colors(cc$palette), length(miss_d))
        names(extra) <- miss_d
        day_cols <- c(day_cols, extra)
    }
    day_cols <- day_cols[day_lv]

    n_day <- length(day_lv)
    n_tr <- length(tr_lv)
    # Divider x positions between treatment blocks (continuous scale).
    div_x <- if (n_tr > 1L) {
        (seq_len(n_tr - 1L) * n_day) + 0.5
    } else {
        numeric()
    }
    # Treatment label centres (in discrete x integer units).
    tr_lab <- data.frame(
        treatment = tr_lv,
        x = (seq_len(n_tr) - 1L) * n_day + (n_day + 1) / 2,
        stringsAsFactors = FALSE
    )

    p <- ggplot(df, aes(x = .x, y = .data[[y_col]])) +
        geom_boxplot(
            aes(fill = .data[[day_col]]),
            width = cc$box_width,
            outlier.shape = NA,
            colour = "grey20",
            alpha = cc$box_alpha,
            linewidth = 0.45
        ) +
        geom_jitter(
            aes(colour = .data[[day_col]]),
            width = cc$jitter_width,
            height = 0,
            size = cc$point_size,
            alpha = cc$point_alpha,
            stroke = 0
        ) +
        scale_fill_manual(values = day_cols, guide = "none") +
        scale_colour_manual(values = day_cols, guide = "none") +
        scale_x_discrete(
            breaks = x_lv,
            labels = rep(day_lv, times = n_tr),
            drop = FALSE
        ) +
        labs(
            title = if (is.null(config$labels$title) ||
                !nzchar(as.character(config$labels$title))) {
                NULL
            } else {
                config$labels$title
            },
            x = NULL,
            y = config$labels$y
        ) +
        theme_classic(base_size = cc$base_size) +
        theme(
            axis.text = element_text(colour = "black"),
            axis.title = element_text(colour = "black"),
            axis.ticks = element_line(colour = "black"),
            panel.background = element_rect(fill = "white", colour = NA),
            plot.background = element_rect(fill = "white", colour = NA),
            panel.grid = element_blank(),
            plot.title = element_text(
                hjust = 0.5, face = "bold",
                margin = margin(b = 6)
            ),
            plot.margin = margin(10, 12, 28, 10)
        )

    if (length(div_x)) {
        p <- p +
            geom_vline(
                xintercept = div_x,
                linetype = cc$divider_linetype,
                colour = cc$divider_colour,
                linewidth = 0.45
            )
    }

    if (identical(cc$y_trans, "log10")) {
        p <- p + scale_y_log10(
            breaks = scales::trans_breaks("log10", function(x) 10^x),
            labels = scales::trans_format("log10", scales::math_format(10^.x))
        )
        y_for_br <- log10(df[[y_col]])
        y_max <- max(y_for_br[is.finite(y_for_br)])
        to_y <- function(v) 10^v
        y_floor <- 10^floor(min(y_for_br[is.finite(y_for_br)]))
    } else if (identical(cc$y_trans, "identity")) {
        y_for_br <- df[[y_col]]
        y_max <- max(y_for_br[is.finite(y_for_br)])
        to_y <- function(v) v
        y_floor <- min(y_for_br[is.finite(y_for_br)])
    } else {
        stop("cohort$y_trans must be log10 or identity", call. = FALSE)
    }

    br_path <- sidecar_path(cc$brackets_file)
    br_ymax <- y_max
    if (!is.null(br_path)) {
        br <- as.data.frame(read_table_auto(br_path), stringsAsFactors = FALSE)
        require_columns(br, list("group1", "group2", "label"))
        br$x1 <- match(as.character(br$group1), x_lv)
        br$x2 <- match(as.character(br$group2), x_lv)
        if (anyNA(br$x1) || anyNA(br$x2)) {
            stop(
                "brackets.tsv group1/group2 must match treatment|day keys.",
                call. = FALSE
            )
        }
        br$span <- abs(br$x2 - br$x1)
        br <- br[order(br$span, br$x1), , drop = FALSE]
        br$y <- y_max + cc$bracket_pad +
            (seq_len(nrow(br)) - 1L) * cc$bracket_step
        tk <- cc$bracket_tick
        seg <- rbind(
            data.frame(x = br$x1, xend = br$x2, y = br$y, yend = br$y),
            data.frame(x = br$x1, xend = br$x1, y = br$y - tk, yend = br$y),
            data.frame(x = br$x2, xend = br$x2, y = br$y - tk, yend = br$y)
        )
        br_ymax <- max(br$y) + cc$bracket_step * 0.55
        p <- p +
            geom_segment(
                data = seg,
                aes(x = x, xend = xend, y = to_y(y), yend = to_y(yend)),
                inherit.aes = FALSE, linewidth = 0.45, colour = "grey20"
            ) +
            geom_text(
                data = br,
                aes(
                    x = (x1 + x2) / 2, y = to_y(y + tk * 0.9),
                    label = label
                ),
                inherit.aes = FALSE, size = cc$label_size, colour = "grey10",
                vjust = 0, fontface = "bold"
            )
    }

    # Underlined treatment labels below the day ticks.
    if (identical(cc$y_trans, "log10")) {
        tr_lab$y <- 10^(min(y_for_br[is.finite(y_for_br)]) - 0.55)
    } else {
        tr_lab$y <- y_floor - 0.08 * diff(range(y_for_br, finite = TRUE))
    }
    # Underline spans the three day positions of each treatment.
    und <- data.frame(
        x = (seq_len(n_tr) - 1L) * n_day + 1,
        xend = seq_len(n_tr) * n_day,
        y = tr_lab$y,
        yend = tr_lab$y,
        stringsAsFactors = FALSE
    )
    p <- p +
        geom_segment(
            data = und,
            aes(x = x, xend = xend, y = y, yend = yend),
            inherit.aes = FALSE,
            linewidth = 0.55, colour = "grey20"
        ) +
        geom_text(
            data = tr_lab,
            aes(x = x, y = y, label = treatment),
            inherit.aes = FALSE,
            vjust = 1.6,
            size = cc$treatment_label_size,
            colour = "grey10"
        ) +
        coord_cartesian(
            ylim = c(NA, to_y(br_ymax)),
            clip = "off"
        )

    if (!is.null(config$labels$title) && nzchar(as.character(config$labels$title))) {
        p <- p + theme(
            plot.title = element_text(
                hjust = 0.5, face = "bold",
                margin = margin(b = 10)
            ),
            axis.line.x = element_line(colour = "black", linewidth = 0.45)
        )
    }

    save_ggplot(p, io$output, width = cc$width, height = cc$height)
    quit(save = "no")
}

# -----------------------------------------------------------------------------
# PLOT — beeswarm / notch / basic
# -----------------------------------------------------------------------------
require_columns(df, config$columns[c("x", "y")])
x_col <- config$columns$x
y_col <- config$columns$y
df[[y_col]] <- as.numeric(df[[y_col]])
lv <- unique(as.character(df[[x_col]]))
if (!is.null(config$group_order)) {
    lv <- c(
        intersect(as.character(config$group_order), lv),
        setdiff(lv, as.character(config$group_order))
    )
}
df[[x_col]] <- factor(df[[x_col]], levels = lv)
df <- df[is.finite(df[[y_col]]) & !is.na(df[[x_col]]), , drop = FALSE]
if (!nrow(df)) stop("No finite rows to plot.", call. = FALSE)

if (identical(config$style, "beeswarm")) {
    load_packages(c("ggbeeswarm"))
    bc <- config$beeswarm
    n_g <- nlevels(df[[x_col]])
    point_cols <- expand_palette(palette_colors(bc$palette), n_g)
    names(point_cols) <- levels(df[[x_col]])
    median_col <- pal_ref(bc$median)

    x_labels <- levels(df[[x_col]])
    lab_path <- sidecar_path(bc$labels_file)
    if (!is.null(lab_path)) {
        lab <- as.data.frame(read_table_auto(lab_path), stringsAsFactors = FALSE)
        require_columns(lab, list(x_col, "label"))
        lab_map <- stats::setNames(as.character(lab$label), as.character(lab[[x_col]]))
        x_labels <- vapply(levels(df[[x_col]]), function(g) {
            if (!is.null(lab_map[[g]]) && nzchar(lab_map[[g]])) {
                gsub("|", "\n", lab_map[[g]], fixed = TRUE)
            } else {
                g
            }
        }, character(1), USE.NAMES = FALSE)
    }

    p <- ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]])) +
        geom_boxplot(
            width = 0.55, outlier.shape = NA, fill = NA,
            colour = bc$box_colour, linewidth = bc$box_linewidth,
            median.colour = median_col, median.linewidth = 1.1
        ) +
        ggbeeswarm::geom_quasirandom(
            aes(colour = .data[[x_col]]),
            width = bc$point_width, size = bc$point_size,
            alpha = bc$point_alpha, stroke = 0
        ) +
        scale_colour_manual(values = point_cols, guide = "none") +
        scale_x_discrete(labels = x_labels) +
        labs(
            title = if (is.null(config$labels$title) ||
                !nzchar(as.character(config$labels$title))) {
                NULL
            } else {
                config$labels$title
            },
            x = config$labels$x,
            y = config$labels$y
        ) +
        theme_classic(base_size = bc$base_size) +
        theme(
            axis.text = element_text(colour = "black"),
            axis.title = element_text(colour = "black"),
            axis.text.x = element_text(size = 9, lineheight = 0.95),
            axis.ticks = element_line(colour = "black"),
            panel.background = element_rect(fill = "white", colour = NA),
            plot.background = element_rect(fill = "white", colour = NA),
            panel.grid.major = if (isTRUE(bc$grid)) {
                element_line(colour = "grey88", linewidth = 0.3)
            } else {
                element_blank()
            },
            panel.grid.minor = element_blank()
        )

    if (identical(bc$y_trans, "log10")) {
        p <- p + scale_y_log10(
            breaks = scales::trans_breaks("log10", function(x) 10^x),
            labels = scales::trans_format("log10", scales::math_format(10^.x))
        )
        y_for_br <- log10(df[[y_col]])
        y_max <- max(y_for_br[is.finite(y_for_br)])
        to_y <- function(v) 10^v
    } else if (identical(bc$y_trans, "identity")) {
        y_for_br <- df[[y_col]]
        y_max <- max(y_for_br[is.finite(y_for_br)])
        to_y <- function(v) v
    } else {
        stop("beeswarm$y_trans must be log10 or identity", call. = FALSE)
    }

    br_path <- sidecar_path(bc$brackets_file)
    if (!is.null(br_path)) {
        br <- as.data.frame(read_table_auto(br_path), stringsAsFactors = FALSE)
        require_columns(br, list("group1", "group2", "label"))
        lv <- levels(df[[x_col]])
        br$x1 <- match(as.character(br$group1), lv)
        br$x2 <- match(as.character(br$group2), lv)
        if (anyNA(br$x1) || anyNA(br$x2)) {
            stop("brackets.tsv groups must match the x column levels.", call. = FALSE)
        }
        # Stack brackets by span so shorter comparisons sit lower.
        br$span <- abs(br$x2 - br$x1)
        br <- br[order(br$span, br$x1), , drop = FALSE]
        br$y <- y_max + bc$bracket_pad +
            (seq_len(nrow(br)) - 1L) * bc$bracket_step
        tk <- bc$bracket_tick
        seg <- rbind(
            data.frame(x = br$x1, xend = br$x2, y = br$y, yend = br$y),
            data.frame(x = br$x1, xend = br$x1, y = br$y - tk, yend = br$y),
            data.frame(x = br$x2, xend = br$x2, y = br$y - tk, yend = br$y)
        )
        p <- p +
            geom_segment(
                data = seg,
                aes(x = x, xend = xend, y = to_y(y), yend = to_y(yend)),
                inherit.aes = FALSE, linewidth = 0.45, colour = "grey20"
            ) +
            geom_text(
                data = br,
                aes(
                    x = (x1 + x2) / 2, y = to_y(y + tk * 0.9),
                    label = label
                ),
                inherit.aes = FALSE, size = bc$label_size, colour = "grey10",
                vjust = 0
            ) +
            coord_cartesian(
                ylim = c(NA, to_y(max(br$y) + bc$bracket_step * 0.55)),
                clip = "off"
            )
    }

    # Empty x title when labels already carry the construct text.
    if (is.null(config$labels$x) || !nzchar(as.character(config$labels$x))) {
        p <- p + theme(axis.title.x = element_blank())
    }
    if (is.null(config$labels$title) || !nzchar(as.character(config$labels$title))) {
        p <- p + theme(plot.title = element_blank())
    }

    save_ggplot(p, io$output, width = bc$width, height = bc$height)
    quit(save = "no")
}

if (identical(config$style, "notch")) {
    nc <- config$notch
    cols <- expand_palette(
        palette_colors(nc$palette), nlevels(df[[x_col]])
    )
    names(cols) <- levels(df[[x_col]])
    horiz <- identical(nc$orientation, "h")
    if (horiz) {
        # Top of the plot = first group_order level.
        df[[x_col]] <- factor(df[[x_col]], levels = rev(levels(df[[x_col]])))
        cols <- cols[levels(df[[x_col]])]
        p <- ggplot(df, aes(
            x = .data[[y_col]],
            y = .data[[x_col]],
            colour = .data[[x_col]]
        ))
        jitter_w <- 0
        jitter_h <- nc$jitter_width
    } else {
        p <- ggplot(df, aes(
            x = .data[[x_col]],
            y = .data[[y_col]],
            colour = .data[[x_col]]
        ))
        jitter_w <- nc$jitter_width
        jitter_h <- 0
    }
    p <- p +
        geom_boxplot(
            notch = isTRUE(nc$notch),
            notchwidth = nc$notchwidth,
            fill = nc$fill,
            alpha = nc$fill_alpha,
            outlier.shape = nc$outlier.shape,
            linewidth = nc$box_linewidth,
            show.legend = FALSE
        ) +
        geom_jitter(
            width = jitter_w,
            height = jitter_h,
            size = nc$point_size,
            alpha = nc$point_alpha,
            show.legend = FALSE
        ) +
        scale_colour_manual(values = cols, guide = "none") +
        labs(
            title = if (is.null(config$labels$title) ||
                !nzchar(as.character(config$labels$title))) {
                NULL
            } else {
                config$labels$title
            },
            x = config$labels$x,
            y = config$labels$y
        ) +
        theme_classic(base_size = nc$base_size) +
        theme(
            axis.text = element_text(colour = "black"),
            axis.ticks = element_line(colour = "black"),
            panel.background = element_rect(fill = "white", colour = NA),
            plot.background = element_rect(fill = "white", colour = NA),
            panel.grid = element_blank(),
            plot.title = element_text(hjust = 0.5, face = "bold")
        )
    if (is.null(config$labels$title) || !nzchar(as.character(config$labels$title))) {
        p <- p + theme(plot.title = element_blank())
    }
    if (is.null(config$labels$y) || !nzchar(as.character(config$labels$y))) {
        p <- p + theme(axis.title.y = element_blank())
    }
    save_ggplot(p, io$output, width = nc$width, height = nc$height)
    quit(save = "no")
}

fill_values <- c("steelblue", "orange", "skyblue", "red", "green")

p <- ggplot(df, aes(
    x = .data[[x_col]],
    y = .data[[y_col]]
)) +
    geom_boxplot(aes(fill = .data[[x_col]]), alpha = 0.7) +
    geom_jitter(
        colour = "grey50",
        size = 2,
        alpha = 0.5,
        width = 0.2,
        height = 0
    ) +
    scale_fill_manual(values = fill_values) +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y
    ) +
    theme_prism(axis_text_angle = 45) +
    theme(plot.background = element_blank())

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output)
