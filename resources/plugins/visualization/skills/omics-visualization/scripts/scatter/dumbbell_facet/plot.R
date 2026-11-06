#!/usr/bin/env Rscript

# Template-ID: scatter-dumbbell-facet
#
# Purpose:
#   Draw per-sample vertical dumbbells faceted by a clinical response
#   (or other panel factor): connected series points with distinct
#   shapes/colours, panel fills, and a signed / mirrored numeric axis
#   (preview.png; e.g. maxMAF by BOR).
#
# Inputs:
#   One row per sample-by-series. Default example:
#     - sample, response, series, value, line_class
#   value is signed for mirrored display: baseline > 0 (upper half),
#   on-therapy MAF < 0 (lower half), ctDNA clearance (mR) at 0.
#   Tick labels use abs(value).
#
# Output:
#   A PDF, PNG, or SVG faceted sample dumbbell plot.
#
# Dependencies:
#   ggplot2, readr, ggprism, scales, ggnewscale
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   Edit CONFIG for columns, facet/series order, shapes/colours,
#   panel fills, y breaks, and labels. Extra samples or series are
#   extra rows/levels, not a new id. A plain two-point dumbbell stays
#   scatter-dumbbell; multi-group horizontal stays
#   scatter-dumbbell-groups.
#
# Scientific assumptions:
#   value, series, and response are supplied. This script does not
#   compute MAF, call response, clearance, or test differences.
#   Sign of value is a display convention (upper = baseline,
#   lower = on-therapy, 0 = mR clearance mark).

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
        x = "sample",
        y = "value",
        facet = "response",
        series = "series",
        line_class = "line_class"
    ),
    facet_order = c("CR", "PR", "SD", "PD"),
    series_order = c("Baseline", "mR", "C2D1", "C3D1"),
    # Panel fill by facet level (unlisted levels use grey95).
    panel_fill = c(
        CR = "#D6EAF8",
        PR = "#D6EAF8",
        SD = "#F5D6D6",
        PD = "#F5D6D6"
    ),
    # Strip (header) fill by facet level.
    strip_fill = c(
        CR = "#5B9BD5",
        PR = "#5B9BD5",
        SD = "#C0504D",
        PD = "#C0504D"
    ),
    series_colours = c(
        Baseline = "#2CA02C",
        mR = "#E36C2C",
        C2D1 = "#7E57C2",
        C3D1 = "#F0C419"
    ),
    series_shapes = c(
        Baseline = 16,
        mR = 17,
        C2D1 = 16,
        C3D1 = 16
    ),
    line_colours = c(
        benefit = "#5B9BD5",
        nonbenefit = "#C4783A"
    ),
    y_breaks = c(-60, -20, -10, -5, -1, 0, 1, 5, 10, 20, 60),
    y_limits = c(-60, 60),
    y_sigma = 1.2,
    point_size = 2.8,
    line_width = 0.55,
    labels = list(
        title = NULL,
        x = NULL,
        y = "maxMAF (%)"
    ),
    size = list(width = 11.5, height = 5.6)
)

load_packages(c("ggplot2", "readr", "ggprism", "scales", "ggnewscale"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- as.data.frame(read_table_auto(io$input), stringsAsFactors = FALSE)
require_columns(df, config$columns)
x_col <- config$columns$x
y_col <- config$columns$y
f_col <- config$columns$facet
s_col <- config$columns$series
l_col <- config$columns$line_class

df[[y_col]] <- as.numeric(df[[y_col]])
df <- df[is.finite(df[[y_col]]), , drop = FALSE]
if (!nrow(df)) stop("No finite value rows.", call. = FALSE)

f_lv <- config$facet_order
f_lv <- c(intersect(f_lv, unique(as.character(df[[f_col]]))),
          setdiff(unique(as.character(df[[f_col]])), f_lv))
s_lv <- config$series_order
s_lv <- c(intersect(s_lv, unique(as.character(df[[s_col]]))),
          setdiff(unique(as.character(df[[s_col]])), s_lv))

df[[f_col]] <- factor(df[[f_col]], levels = f_lv)
df[[s_col]] <- factor(df[[s_col]], levels = s_lv)
df[[l_col]] <- factor(as.character(df[[l_col]]))
# Unique x within facet so free_x does not repeat every sample in every panel.
df$.x_id <- interaction(
    df[[f_col]], as.character(df[[x_col]]),
    sep = "\r", drop = TRUE, lex.order = FALSE
)
x_lab <- setNames(
    as.character(df[[x_col]]),
    as.character(df$.x_id)
)
x_lab <- x_lab[!duplicated(names(x_lab))]

seg <- do.call(rbind, lapply(
    split(df, df$.x_id, drop = TRUE),
    function(d) {
        out <- data.frame(
            .x_id = d$.x_id[[1]],
            line = d[[l_col]][[1]],
            ymin = min(d[[y_col]]),
            ymax = max(d[[y_col]]),
            stringsAsFactors = FALSE
        )
        out[[f_col]] <- d[[f_col]][[1]]
        out
    }
))
seg[[f_col]] <- factor(seg[[f_col]], levels = levels(df[[f_col]]))
seg$.x_id <- factor(as.character(seg$.x_id), levels = levels(df$.x_id))

# Panel background rectangles in data space (one per facet).
panel_df <- data.frame(
    xmin = rep(-Inf, length(f_lv)),
    xmax = rep(Inf, length(f_lv)),
    ymin = rep(-Inf, length(f_lv)),
    ymax = rep(Inf, length(f_lv)),
    stringsAsFactors = FALSE
)
panel_df[[f_col]] <- factor(f_lv, levels = f_lv)
panel_df$fill <- ifelse(
    as.character(panel_df[[f_col]]) %in% names(config$panel_fill),
    unname(config$panel_fill[as.character(panel_df[[f_col]])]),
    "#F5F5F5"
)

# Series aesthetics (expand if extra series appear).
ser_cols <- config$series_colours
ser_shp <- config$series_shapes
miss_s <- setdiff(s_lv, names(ser_cols))
if (length(miss_s)) {
    extra <- expand_palette(palette_colors("Qualitative.Bold"), length(miss_s))
    names(extra) <- miss_s
    ser_cols <- c(ser_cols, extra)
}
miss_sh <- setdiff(s_lv, names(ser_shp))
if (length(miss_sh)) {
    ser_shp <- c(ser_shp, setNames(rep(16L, length(miss_sh)), miss_sh))
}
ser_cols <- ser_cols[s_lv]
ser_shp <- ser_shp[s_lv]

line_cols <- config$line_colours
miss_l <- setdiff(levels(df[[l_col]]), names(line_cols))
if (length(miss_l)) {
    extra <- expand_palette(palette_colors("Qualitative.Safe"), length(miss_l))
    names(extra) <- miss_l
    line_cols <- c(line_cols, extra)
}

y_breaks <- config$y_breaks
y_labs <- as.character(abs(y_breaks))

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p <- ggplot() +
    geom_rect(
        data = panel_df,
        aes(
            xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
            fill = .data[[f_col]]
        ),
        colour = NA,
        show.legend = FALSE
    ) +
    scale_fill_manual(
        values = setNames(panel_df$fill, as.character(panel_df[[f_col]])),
        guide = "none"
    ) +
    geom_hline(
        yintercept = 0,
        linetype = "dashed",
        colour = "grey40",
        linewidth = 0.4
    ) +
    geom_segment(
        data = seg,
        aes(x = .x_id, xend = .x_id, y = ymin, yend = ymax, colour = line),
        linewidth = config$line_width,
        show.legend = FALSE
    ) +
    scale_colour_manual(values = line_cols, guide = "none") +
    ggnewscale::new_scale_colour() +
    geom_point(
        data = df,
        aes(
            x = .x_id,
            y = .data[[y_col]],
            colour = .data[[s_col]],
            shape = .data[[s_col]]
        ),
        size = config$point_size,
        stroke = 0.9
    ) +
    facet_grid(
        cols = vars(.data[[f_col]]),
        scales = "free_x",
        space = "free_x"
    ) +
    scale_colour_manual(values = ser_cols, name = NULL) +
    scale_shape_manual(values = ser_shp, name = NULL) +
    scale_y_continuous(
        trans = scales::pseudo_log_trans(base = 10, sigma = config$y_sigma),
        breaks = y_breaks,
        labels = y_labs,
        limits = config$y_limits,
        expand = expansion(mult = c(0.04, 0.06))
    ) +
    scale_x_discrete(
        expand = expansion(add = 0.6),
        labels = function(x) {
            unname(x_lab[x])
        }
    ) +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y
    ) +
    theme_prism(base_size = 11) +
    theme(
        panel.grid = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = rel(0.55)),
        axis.text.y = element_text(size = rel(0.85)),
        strip.background = element_rect(colour = NA),
        strip.text.x = element_text(
            face = "bold", colour = "white", size = rel(1.0),
            hjust = 0.5, vjust = 0.5, margin = margin(t = 5, b = 5)
        ),
        legend.position = "top",
        legend.direction = "horizontal",
        panel.spacing.x = unit(0.15, "lines"),
        plot.background = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = NA, colour = NA),
        plot.margin = margin(6, 8, 4, 6)
    )

# Per-facet strip colours + force vertical centre of strip labels.
gt <- ggplotGrob(p)
strip_idx <- which(grepl("strip-t", gt$layout$name))
centre_strip_text <- function(g) {
    if (inherits(g, "text")) {
        # Slightly below geometric centre so bold caps look optically centred.
        g$y <- unit(0.46, "npc")
        g$vjust <- 0.5
        return(g)
    }
    if (!is.null(g$children) && length(g$children)) {
        for (k in seq_along(g$children)) {
            g$children[[k]] <- centre_strip_text(g$children[[k]])
        }
    }
    if (!is.null(g$grobs) && length(g$grobs)) {
        for (k in seq_along(g$grobs)) {
            g$grobs[[k]] <- centre_strip_text(g$grobs[[k]])
        }
    }
    g
}
for (i in seq_along(strip_idx)) {
    ii <- strip_idx[[i]]
    fac <- as.character(f_lv[[min(i, length(f_lv))]])
    col <- if (fac %in% names(config$strip_fill)) {
        config$strip_fill[[fac]]
    } else {
        "grey50"
    }
    if (length(gt$grobs[[ii]]$grobs)) {
        for (j in seq_along(gt$grobs[[ii]]$grobs)) {
            g <- gt$grobs[[ii]]$grobs[[j]]
            if (inherits(g, "rect") || inherits(g, "roundrect")) {
                gt$grobs[[ii]]$grobs[[j]]$gp$fill <- col
            }
            if (!is.null(g$children)) {
                for (k in seq_along(g$children)) {
                    ch <- g$children[[k]]
                    if (inherits(ch, "rect") || inherits(ch, "roundrect")) {
                        gt$grobs[[ii]]$grobs[[j]]$children[[k]]$gp$fill <- col
                    }
                }
            }
        }
    }
    gt$grobs[[ii]] <- centre_strip_text(gt$grobs[[ii]])
}

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(
    gt,
    io$output,
    width = config$size$width,
    height = config$size$height
)
