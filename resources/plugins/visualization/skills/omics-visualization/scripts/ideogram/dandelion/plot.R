#!/usr/bin/env Rscript

# Template-ID: ideogram-dandelion
#
# Purpose:
#   Draw dense sites as dandelions: neighbouring positions collapse
#   into one stem whose height is cluster size (or mean score), with
#   fan / pie / circle / pin heads and a feature bar.
#
# Inputs:
#   One row per site. Default example:
#     - pos: genomic or local coordinate
#     - score: supplied value (fan/pie expect [0, 1])
#   Optional sidecar beside the input:
#     - features.tsv: start, end, feature
#
# Output:
#   A PDF, PNG, or SVG dandelion figure.
#
# Dependencies:
#   ggplot2, readr, ggprism, patchwork
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (columns, type, maxgaps,
#   height method). pin / fan / pie / circle stay on this id. Sparse
#   one-stem-per-site lollipops are ideogram-lollipop. Whole-chromosome
#   window fill is ideogram-density. Do not fetch TxDb, UCSC, or VCF.
#
# Scientific assumptions:
#   Positions and scores are supplied. Grouping sites closer than
#   maxgaps is a display aggregate, not a peak caller. Stem height
#   n is the count of supplied rows in that group; mean is the mean
#   of supplied scores. This script does not compute methylation, call
#   SNPs, or load a gene model.

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
        pos = "pos",
        score = "score",
        score2 = "score2"
    ),
    sidecars = list(
        features = "features.tsv"
    ),
    feature = list(
        start = "start",
        end = "end",
        name = "feature"
    ),
    type = "pin",
    maxgaps = 1 / 25,
    height_method = "n",
    fig = list(
        width = 8.6,
        height = 4.4,
        panel_width = 7.4,
        panel_height = 2.9
    ),
    labels = list(
        title = NULL,
        x = "Position",
        y = NULL,
        fill = NULL
    )
)

load_packages(c("ggplot2", "readr", "ggprism", "patchwork"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
need <- config$columns[c("pos", "score")]
require_columns(df, need)

pos_col <- config$columns$pos
sc_col <- config$columns$score
s2_col <- config$columns$score2

df[[pos_col]] <- as.numeric(df[[pos_col]])
df[[sc_col]] <- as.numeric(df[[sc_col]])
ok <- is.finite(df[[pos_col]]) & is.finite(df[[sc_col]])
if (any(!ok)) {
    message("Dropped ", sum(!ok), " site row(s) with invalid position or score.")
    df <- df[ok, , drop = FALSE]
}
if (!nrow(df)) {
    stop("No site rows left to plot.", call. = FALSE)
}

head_type <- tolower(config$type)
if (!head_type %in% c("fan", "circle", "pie", "pin")) {
    stop("config$type must be fan, circle, pie, or pin.", call. = FALSE)
}
hm <- tolower(config$height_method)
if (!hm %in% c("n", "mean")) {
    stop("config$height_method must be \"n\" or \"mean\".", call. = FALSE)
}

x_min <- min(df[[pos_col]])
x_max <- max(df[[pos_col]])

side_path <- function(name) {
    if (is.null(name) || !nzchar(name)) {
        return(NA_character_)
    }
    file.path(dirname(io$input), name)
}

features <- NULL
feat_path <- side_path(config$sidecars$features)
if (!is.na(feat_path) && file.exists(feat_path)) {
    features <- read_table_auto(feat_path)
    require_columns(features, config$feature)
    fs <- config$feature$start
    fe <- config$feature$end
    fn <- config$feature$name
    features[[fs]] <- as.numeric(features[[fs]])
    features[[fe]] <- as.numeric(features[[fe]])
    ok_f <- is.finite(features[[fs]]) & is.finite(features[[fe]]) &
        !is.na(features[[fn]])
    if (any(!ok_f)) {
        features <- features[ok_f, , drop = FALSE]
    }
    if (!nrow(features)) {
        features <- NULL
    } else {
        swap_f <- features[[fs]] > features[[fe]]
        if (any(swap_f)) {
            tmp <- features[[fs]][swap_f]
            features[[fs]][swap_f] <- features[[fe]][swap_f]
            features[[fe]][swap_f] <- tmp
        }
        features[[fn]] <- factor(
            as.character(features[[fn]]),
            levels = unique(as.character(features[[fn]]))
        )
        x_min <- min(x_min, min(features[[fs]]))
        x_max <- max(x_max, max(features[[fe]]))
    }
}

span <- x_max - x_min
if (!is.finite(span) || span <= 0) {
    span <- 1
}
mg <- config$maxgaps
if (is.null(mg) || !is.finite(mg) || mg <= 0) {
    mg <- 1 / 50
}
max_gap <- if (mg <= 1) span * mg else mg

o <- order(df[[pos_col]])
pos_o <- df[[pos_col]][o]
gps_o <- cumsum(c(TRUE, diff(pos_o) > max_gap))
df$gps <- NA_integer_
df$gps[o] <- gps_o

score_raw <- df[[sc_col]]
if (identical(head_type, "fan") && any(score_raw > 1, na.rm = TRUE)) {
    message("Fan heads expect score in [0, 1]; dividing by the maximum as a display scale.")
    df[[sc_col]] <- score_raw / max(score_raw, na.rm = TRUE)
}
if (is.null(s2_col) || !nzchar(s2_col) || !(s2_col %in% names(df))) {
    df$score2 <- pmax(0, 1 - df[[sc_col]])
    s2_col <- "score2"
} else {
    df[[s2_col]] <- as.numeric(df[[s2_col]])
    miss2 <- !is.finite(df[[s2_col]])
    df[[s2_col]][miss2] <- pmax(0, 1 - df[[sc_col]][miss2])
}

stems <- do.call(rbind, lapply(split(df, df$gps), function(g) {
    data.frame(
        gps = g$gps[[1]],
        n = nrow(g),
        mean_score = mean(g[[sc_col]]),
        mid = mean(g[[pos_col]]),
        lo = min(g[[pos_col]]),
        hi = max(g[[pos_col]]),
        stringsAsFactors = FALSE
    )
}))
stems$height <- if (identical(hm, "mean")) stems$mean_score else stems$n
max_h <- max(stems$height)
if (!is.finite(max_h) || max_h <= 0) {
    stop("Cluster heights are not finite.", call. = FALSE)
}

y_lab <- config$labels$y
if (is.null(y_lab)) {
    y_lab <- if (identical(hm, "mean")) "Mean score" else "Sites per cluster"
}

df$mid <- stems$mid[match(df$gps, stems$gps)]
df$height <- stems$height[match(df$gps, stems$gps)]
df$n_cl <- stems$n[match(df$gps, stems$gps)]
df$lo <- stems$lo[match(df$gps, stems$gps)]
df$hi <- stems$hi[match(df$gps, stems$gps)]

x_span <- span
y_span <- max_h * 1.28
# Puff size follows cluster n (a round seed head), not genomic width.
# Wide clusters were drawing triangular trunks; tight sites + this
# radius keep a thin stalk and a circular crown.
df$r_y <- 1.55 + 0.82 * sqrt(df$n_cl)
df$r_x <- df$r_y * (x_span / y_span) *
    (config$fig$panel_height / config$fig$panel_width)
y_top <- max(df$height + df$r_y) * 1.18
y_span <- y_top
df$r_x <- df$r_y * (x_span / y_span) *
    (config$fig$panel_height / config$fig$panel_width)

pin_h_mm <- 3.8
pin_full_y <- (pin_h_mm / 25.4) / config$fig$panel_height * y_span
pin_full_x <- pin_full_y * (x_span / y_span) *
    (config$fig$panel_height / config$fig$panel_width)
pin_rx <- 0.34 * pin_full_x
pin_ry <- 0.34 * pin_full_y
pin_lift <- 0.66 * pin_full_y

puff <- vector("list", nrow(df))
for (g in unique(df$gps)) {
    idx <- which(df$gps == g)
    n <- length(idx)
    mid <- df$mid[idx[[1]]]
    h <- df$height[idx[[1]]]
    rx <- df$r_x[idx[[1]]]
    ry <- df$r_y[idx[[1]]]
    if (n == 1L) {
        i <- idx[[1]]
        puff[[i]] <- data.frame(
            x0 = mid,
            y0 = 0,
            x1 = mid,
            y1 = h,
            x2 = mid,
            y2 = h,
            rot = 0,
            stringsAsFactors = FALSE
        )
    } else {
        # Semicircle crown (trackViewer pin / fan heads sit on this arc).
        theta <- seq(-0.96, 0.96, length.out = n) * (pi / 2)
        for (j in seq_along(idx)) {
            i <- idx[[j]]
            puff[[i]] <- data.frame(
                x0 = mid,
                y0 = 0,
                x1 = mid,
                y1 = h,
                x2 = mid + rx * sin(theta[[j]]),
                y2 = h + ry * cos(theta[[j]]),
                rot = theta[[j]],
                stringsAsFactors = FALSE
            )
        }
    }
}
puff <- do.call(rbind, puff)
df$x0 <- puff$x0
df$y0 <- puff$y0
df$x1 <- puff$x1
df$y1 <- puff$y1
df$x2 <- puff$x2
df$y2 <- puff$y2
df$rot <- puff$rot
df$percent <- pmin(pmax(df[[sc_col]], 0), 1)

bold <- palette_colors("Qualitative.Bold", n = 4)
fan_fill <- bold[[2]]
pin_body <- bold[[2]]
pin_inner <- bold[[4]]
pie_cols <- bold[c(2, 4)]
names(pie_cols) <- c("score", "other")
feat_cols <- NULL
if (!is.null(features)) {
    fn <- config$feature$name
    n_f <- nlevels(features[[fn]])
    if (n_f == 2L) {
        feat_cols <- palette_colors("Qualitative.Paired", n = 5)[c(1, 5)]
    } else {
        idx <- seq(1L, by = 2L, length.out = n_f)
        feat_cols <- palette_colors("Qualitative.Paired", n = max(idx))[idx]
    }
    names(feat_cols) <- levels(features[[fn]])
}

stem_df <- unique(df[, c("x0", "y0", "x1", "y1")])
y_top <- max(df$y2) + 1.15 * pin_full_y
pad <- 0.04 * span
x_limits <- c(
    min(x_min, min(df$x2) - pin_full_x) - pad,
    max(x_max, max(df$x2) + pin_full_x) + pad
)

pin_poly <- NULL
pin_dot <- NULL
if (identical(head_type, "pin")) {
    pin_rows <- vector("list", nrow(df))
        th <- seq(pi * 1.18, -pi * 0.18, length.out = 30)
    for (i in seq_len(nrow(df))) {
        cx <- df$x2[[i]]
        cy <- df$y2[[i]] + pin_lift
        pin_rows[[i]] <- data.frame(
            x = c(df$x2[[i]], cx + pin_rx * cos(th)),
            y = c(df$y2[[i]], cy + pin_ry * sin(th)),
            grp = i
        )
    }
    pin_poly <- do.call(rbind, pin_rows)
    pin_dot <- data.frame(
        x = df$x2,
        y = df$y2 + pin_lift
    )
}

fan_poly <- NULL
if (identical(head_type, "fan")) {
    fan_rows <- vector("list", nrow(df))
    n_edge <- 18L
    for (i in seq_len(nrow(df))) {
        pct <- df$percent[[i]]
        if (pct <= 0) {
            next
        }
        rx <- 0.22 * df$r_x[[i]] * (0.55 + 0.45 * pct)
        ry <- 0.22 * df$r_y[[i]] * (0.55 + 0.45 * pct)
        th <- seq(-pct / 2, pct / 2, length.out = n_edge) * pi +
            pi / 2 - df$rot[[i]]
        fan_rows[[i]] <- data.frame(
            x = c(df$x2[[i]], df$x2[[i]] + rx * cos(th)),
            y = c(df$y2[[i]], df$y2[[i]] + ry * sin(th)),
            grp = i
        )
    }
    fan_poly <- do.call(rbind, fan_rows)
}

pie_poly <- NULL
if (identical(head_type, "pie")) {
    pie_rows <- vector("list", nrow(df))
    for (i in seq_len(nrow(df))) {
        tot <- df[[sc_col]][[i]] + df[[s2_col]][[i]]
        if (!is.finite(tot) || tot <= 0) {
            next
        }
        f1 <- df[[sc_col]][[i]] / tot
        rx <- 0.22 * df$r_x[[i]]
        ry <- 0.22 * df$r_y[[i]]
        th1 <- seq(0, f1, length.out = 18) * 2 * pi + pi / 2
        th2 <- seq(f1, 1, length.out = 18) * 2 * pi + pi / 2
        pie_rows[[i]] <- rbind(
            data.frame(
                x = c(df$x2[[i]], df$x2[[i]] + rx * cos(th1), df$x2[[i]]),
                y = c(df$y2[[i]], df$y2[[i]] + ry * sin(th1), df$y2[[i]]),
                grp = paste0(i, "a"),
                slice = "score"
            ),
            data.frame(
                x = c(df$x2[[i]], df$x2[[i]] + rx * cos(th2), df$x2[[i]]),
                y = c(df$y2[[i]], df$y2[[i]] + ry * sin(th2), df$y2[[i]]),
                grp = paste0(i, "b"),
                slice = "other"
            )
        )
    }
    pie_poly <- do.call(rbind, pie_rows)
}

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
p_main <- ggplot2::ggplot() +
    ggplot2::geom_segment(
        data = stem_df,
        ggplot2::aes(x = x0, xend = x1, y = y0, yend = y1),
        colour = "grey40",
        linewidth = 0.55,
        lineend = "butt"
    ) +
    ggplot2::geom_segment(
        data = df,
        ggplot2::aes(x = x1, xend = x2, y = y1, yend = y2),
        colour = "grey42",
        linewidth = 0.42,
        lineend = "butt"
    )

if (identical(head_type, "fan") && !is.null(fan_poly)) {
    p_main <- p_main +
        ggplot2::geom_polygon(
            data = fan_poly,
            ggplot2::aes(x = x, y = y, group = grp),
            fill = fan_fill,
            colour = "grey25",
            linewidth = 0.15
        )
} else if (identical(head_type, "pie") && !is.null(pie_poly)) {
    p_main <- p_main +
        ggplot2::geom_polygon(
            data = pie_poly,
            ggplot2::aes(x = x, y = y, group = grp, fill = slice),
            colour = "grey25",
            linewidth = 0.12
        ) +
        ggplot2::scale_fill_manual(
            values = pie_cols,
            name = config$labels$fill,
            labels = c(score = "Score", other = "1 \u2212 score")
        )
} else if (identical(head_type, "pin") && !is.null(pin_poly)) {
    p_main <- p_main +
        ggplot2::geom_polygon(
            data = pin_poly,
            ggplot2::aes(x = x, y = y, group = grp),
            fill = pin_body,
            colour = pin_body,
            linewidth = 0.12
        ) +
        ggplot2::geom_point(
            data = pin_dot,
            ggplot2::aes(x = x, y = y),
            shape = 21,
            fill = pin_inner,
            colour = pin_body,
            size = 1.15,
            stroke = 0.25
        )
} else {
    p_main <- p_main +
        ggplot2::geom_point(
            data = df,
            ggplot2::aes(x = x2, y = y2),
            shape = 21,
            fill = fan_fill,
            colour = "grey25",
            size = 2.4,
            stroke = 0.3
        )
}

p_main <- p_main +
    ggplot2::scale_x_continuous(limits = x_limits, expand = c(0, 0)) +
    ggplot2::coord_cartesian(ylim = c(0, y_top), clip = "off") +
    ggplot2::labs(
        title = config$labels$title,
        x = NULL,
        y = y_lab
    ) +
    ggprism::theme_prism() +
    ggplot2::theme(
        plot.background = ggplot2::element_rect(fill = NA, colour = NA),
        panel.background = ggplot2::element_rect(fill = NA, colour = NA),
        legend.background = ggplot2::element_blank(),
        legend.key = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        axis.line.x = ggplot2::element_blank(),
        plot.margin = ggplot2::margin(6, 8, 2, 8)
    )

backbone <- data.frame(start = x_min, end = x_max)
p_feat <- ggplot2::ggplot() +
    ggplot2::geom_rect(
        data = backbone,
        ggplot2::aes(xmin = start, xmax = end, ymin = 0.38, ymax = 0.62),
        fill = "grey78",
        colour = NA
    )
if (!is.null(features)) {
    fs <- config$feature$start
    fe <- config$feature$end
    fn <- config$feature$name
    features$mid <- (features[[fs]] + features[[fe]]) / 2
    p_feat <- p_feat +
        ggplot2::geom_rect(
            data = features,
            ggplot2::aes(
                xmin = .data[[fs]],
                xmax = .data[[fe]],
                ymin = 0.12,
                ymax = 0.88,
                fill = .data[[fn]]
            ),
            colour = "grey20",
            linewidth = 0.25
        ) +
        ggplot2::geom_text(
            data = features,
            ggplot2::aes(x = mid, y = 0.5, label = .data[[fn]]),
            size = 4.2,
            colour = "grey10",
            fontface = "bold"
        ) +
        ggplot2::scale_fill_manual(values = feat_cols, guide = "none")
}
p_feat <- p_feat +
    ggplot2::scale_x_continuous(limits = x_limits, expand = c(0, 0)) +
    ggplot2::coord_cartesian(ylim = c(0, 1), clip = "off") +
    ggplot2::labs(x = config$labels$x, y = NULL) +
    ggprism::theme_prism() +
    ggplot2::theme(
        plot.background = ggplot2::element_rect(fill = NA, colour = NA),
        panel.background = ggplot2::element_rect(fill = NA, colour = NA),
        legend.position = "none",
        axis.line.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        plot.margin = ggplot2::margin(0, 8, 6, 8)
    )

p <- p_main / p_feat + patchwork::plot_layout(heights = c(3.55, 0.72))
g <- patchwork::patchworkGrob(p)
bg <- which(g$layout$name == "background")
if (length(bg)) {
    for (i in bg) {
        g$grobs[[i]]$gp$fill <- NA
        g$grobs[[i]]$gp$col <- NA
    }
}

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(g, io$output, width = config$fig$width, height = config$fig$height)
