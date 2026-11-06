#!/usr/bin/env Rscript

# Template-ID: line-survival
#
# Purpose:
#   Draw a Kaplan-Meier survival curve by group. style = "basic"
#   (default) is a two-group curve with log-rank / HR annotation.
#   style = "multi" draws multi-group curves with CI ribbons,
#   right-side legend, optional overall P, and an optional supplied
#   pairwise significance table inset in open panel space
#   (preview_multi.png). Risk table is off by default.
#
# Inputs:
#   One row per subject. Default example:
#     - OS_time, OS, group
#   Multi example (example_multi.tsv):
#     - time, event, group  (Low / Intermediate / High / Unknown)
#   Optional sidecar (pairwise_multi.tsv):
#     - comparison, p   (pairwise labels are supplied; not computed)
#
# Output:
#   A PDF, PNG, or SVG survival curve.
#
# Dependencies:
#   ggplot2, readr, survival, survminer, ggprism, grid
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   Edit CONFIG for columns, group_order, palette, risk-table breaks,
#   legend, pval_label, and multi$pairwise_file. Extra groups are
#   extra strata, not a new id. Do not add an id for multi-group
#   alone. Patient treatment lanes stay bar-swimmer.
#
# Scientific assumptions:
#   Time and event columns define a right-censored survival object.
#   Rows with missing time, event, or group are dropped. The KM fit
#   is a display of the supplied follow-up. When pval_label is set,
#   that text is used; otherwise a log-rank p is shown. Pairwise
#   comparison text in multi$pairwise_file is display-only.
#   basic-style HR is for the second group relative to the first.

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
    # "basic": classic two-group preview.
    # "multi": multi-strata KM + CI + risk table + right legend.
    style = "basic",
    columns = list(
        time = "OS_time",
        event = "OS",
        group = "group"
    ),
    group_order = NULL,
    labels = list(
        title = "TCGA-BRCA",
        x = "Overall survival (Days)",
        y = "Survival probability",
        legend = "Molecular risk"
    ),
    # NULL = compute log-rank; character = use as annotation text.
    pval_label = NULL,
    palettes = list(
        group = "Qualitative.Bold"
    ),
    multi = list(
        conf_int = TRUE,
        conf_alpha = 0.18,
        risk_table = FALSE,
        risk_table_height = 0.30,
        break_time_by = 6,
        xlim = c(0, 72),
        censor_size = 2.8,
        legend_position = "right",
        # Supplied pairwise table (comparison, p). NULL = omit inset.
        pairwise_file = NULL,
        pairwise_title = "Pairwise log-rank",
        # Keep only rows with p < alpha (NULL = keep all supplied rows).
        pairwise_alpha = 0.05,
        # Overall P (left blank) and pairwise table (right blank).
        pval_pos = c(x = 12, y = 0.06),
        pairwise_box = c(xmin = 40, xmax = 71.5, ymin = 0.00, ymax = 0.25),
        size = list(width = 9.0, height = 6.2)
    ),
    base_size = 12,
    size = list(width = 8, height = 6)
)

load_packages(c("ggplot2", "readr", "survival", "survminer", "ggprism", "grid"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- as.data.frame(read_table_auto(io$input), stringsAsFactors = FALSE)
require_columns(df, config$columns)
time_col <- config$columns$time
event_col <- config$columns$event
group_col <- config$columns$group

df <- df[
    !is.na(df[[time_col]]) &
        !is.na(df[[event_col]]) &
        !is.na(df[[group_col]]),
    ,
    drop = FALSE
]
df[[time_col]] <- as.numeric(df[[time_col]])
df[[event_col]] <- as.numeric(df[[event_col]])
df <- df[is.finite(df[[time_col]]) & is.finite(df[[event_col]]), , drop = FALSE]
if (!nrow(df)) stop("No finite survival rows.", call. = FALSE)

lv <- unique(as.character(df[[group_col]]))
if (!is.null(config$group_order)) {
    lv <- c(
        intersect(as.character(config$group_order), lv),
        setdiff(lv, as.character(config$group_order))
    )
}
df[[group_col]] <- factor(df[[group_col]], levels = lv)

f <- as.formula(sprintf(
    "Surv(%s, %s) ~ `%s`",
    time_col, event_col, group_col
))
fits <- survfit(f, data = df)
fits$call$formula <- f

n_g <- nlevels(df[[group_col]])
cols <- expand_palette(palette_colors(config$palettes$group), n_g)
# Prefer blue / gold / red / grey for Low/Intermediate/High/Unknown when 4+.
bold <- palette_colors("Qualitative.Bold")
safe <- palette_colors("Qualitative.Safe")
if (identical(config$palettes$group, "Qualitative.Bold") && n_g >= 4L) {
    pick <- c(bold[[3]], bold[[4]], bold[[5]], safe[[12]], bold)
    cols <- expand_palette(pick, n_g)
}
names(cols) <- levels(df[[group_col]])

format_p <- function(p) {
    if (!is.finite(p)) return("P = NA")
    if (p < 0.0001) return("P < 0.0001")
    sprintf("P = %.4g", p)
}

pval_text <- config$pval_label
if (is.null(pval_text) || !nzchar(as.character(pval_text))) {
    pval_text <- format_p(surv_pvalue(fits, data = df)$pval)
}

# -----------------------------------------------------------------------------
# MULTI STYLE
# -----------------------------------------------------------------------------
if (identical(config$style, "multi")) {
    mc <- config$multi
    sidecar_path <- function(name) {
        if (is.null(name) || !nzchar(as.character(name)[[1]])) {
            return(NULL)
        }
        name <- as.character(name)[[1]]
        if (file.exists(name)) return(normalizePath(name, mustWork = TRUE))
        cand <- file.path(dirname(io$input), name)
        if (file.exists(cand)) return(normalizePath(cand, mustWork = TRUE))
        file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
        if (length(file_arg)) {
            script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[[1]])))
            cand2 <- file.path(script_dir, name)
            if (file.exists(cand2)) {
                return(normalizePath(cand2, mustWork = TRUE))
            }
        }
        NULL
    }

    th <- theme_prism(base_size = config$base_size) +
        theme(
            panel.grid = element_blank(),
            plot.background = element_rect(fill = "white", colour = NA),
            panel.background = element_rect(fill = "white", colour = NA),
            legend.background = element_blank(),
            legend.title = element_text(face = "bold", size = rel(0.95)),
            legend.text = element_text(size = rel(0.9)),
            plot.title = element_text(hjust = 0.5, face = "bold", size = rel(1.05)),
            axis.title = element_text(size = rel(0.95))
        )

    fit_plot <- ggsurvplot(
        fits,
        data = df,
        conf.int = isTRUE(mc$conf_int),
        conf.int.alpha = mc$conf_alpha,
        censor = TRUE,
        censor.shape = "|",
        censor.size = mc$censor_size,
        pval = FALSE,
        risk.table = isTRUE(mc$risk_table),
        risk.table.col = "strata",
        risk.table.title = "No. at risk",
        risk.table.height = mc$risk_table_height,
        risk.table.y.text = FALSE,
        risk.table.y.text.col = TRUE,
        tables.theme = theme(
            axis.ticks = element_blank(),
            axis.text.x = element_blank(),
            axis.title = element_blank(),
            axis.line = element_blank(),
            panel.grid = element_blank(),
            plot.title = element_text(
                size = rel(0.95), face = "bold", hjust = 0
            ),
            plot.margin = margin(t = 2, r = 8, b = 2, l = 8)
        ),
        legend = mc$legend_position,
        legend.title = config$labels$legend,
        legend.labs = levels(df[[group_col]]),
        palette = unname(cols),
        break.time.by = mc$break_time_by,
        xlim = mc$xlim,
        xlab = config$labels$x,
        ylab = config$labels$y,
        title = config$labels$title,
        ggtheme = th
    )

    # Overall P in the left blank (annotation_custom stays above CI ribbons).
    ppos <- mc$pval_pos
    p_inset <- grid::grobTree(
        grid::rectGrob(
            gp = grid::gpar(fill = "white", col = NA)
        ),
        grid::textGrob(
            pval_text,
            x = unit(0.08, "npc"),
            y = unit(0.5, "npc"),
            just = c("left", "center"),
            gp = grid::gpar(fontsize = 11, fontface = "bold", col = "grey15")
        )
    )
    fit_plot$plot <- fit_plot$plot +
        annotation_custom(
            grob = p_inset,
            xmin = ppos[["x"]],
            xmax = ppos[["x"]] + 14,
            ymin = ppos[["y"]],
            ymax = ppos[["y"]] + 0.08
        )

    pw_path <- sidecar_path(mc$pairwise_file)
    if (!is.null(pw_path)) {
        pw <- as.data.frame(read_table_auto(pw_path), stringsAsFactors = FALSE)
        require_columns(pw, list("comparison", "p"))
        parse_p <- function(x) {
            x <- trimws(as.character(x))
            x <- gsub(",", "", x, fixed = TRUE)
            x <- sub("^[Pp]\\s*[=<>≤≥]\\s*", "", x)
            x <- sub("^<", "", x)
            suppressWarnings(as.numeric(x))
        }
        if (!is.null(mc$pairwise_alpha) && is.finite(mc$pairwise_alpha)) {
            # "<0.0001" parses as 0.0001 (upper bound); keep if < alpha.
            p_num <- parse_p(pw$p)
            keep <- is.finite(p_num) & p_num < mc$pairwise_alpha
            pw <- pw[keep, , drop = FALSE]
        }
        if (!nrow(pw)) {
            stop(
                "No significant pairwise rows after pairwise_alpha filter.",
                call. = FALSE
            )
        }
        # Two columns at fixed x so P header and values share one left edge.
        x_cmp <- unit(0.06, "npc")
        x_p <- unit(0.65, "npc")
        n_row <- nrow(pw)
        y_hdr <- 0.70
        y_top <- 0.55
        y_bot <- 0.14
        ys <- if (n_row == 1L) {
            y_top
        } else {
            y_top - (seq_len(n_row) - 1L) * (y_top - y_bot) / (n_row - 1L)
        }
        grobs <- list(
            grid::rectGrob(
                gp = grid::gpar(fill = "white", col = "grey50", lwd = 0.9)
            ),
            grid::rectGrob(
                x = unit(0.5, "npc"),
                y = unit(y_hdr, "npc"),
                width = unit(0.92, "npc"),
                height = unit(0.14, "npc"),
                gp = grid::gpar(fill = "grey92", col = NA)
            ),
            grid::textGrob(
                as.character(mc$pairwise_title),
                x = x_cmp,
                y = unit(0.88, "npc"),
                just = c("left", "center"),
                gp = grid::gpar(
                    fontsize = 11,
                    fontface = "bold",
                    fontfamily = "Helvetica",
                    col = "black"
                )
            ),
            grid::textGrob(
                "Comparison",
                x = x_cmp,
                y = unit(y_hdr, "npc"),
                just = c("left", "center"),
                gp = grid::gpar(
                    fontsize = 7.5,
                    fontface = "bold",
                    fontfamily = "Helvetica",
                    col = "grey10"
                )
            ),
            grid::textGrob(
                "P",
                x = x_p,
                y = unit(y_hdr, "npc"),
                just = c("left", "center"),
                gp = grid::gpar(
                    fontsize = 7.5,
                    fontface = "bold",
                    fontfamily = "Helvetica",
                    col = "grey10"
                )
            )
        )
        for (i in seq_len(n_row)) {
            yi <- ys[[i]]
            grobs[[length(grobs) + 1L]] <- grid::textGrob(
                as.character(pw$comparison[[i]]),
                x = x_cmp,
                y = unit(yi, "npc"),
                just = c("left", "center"),
                gp = grid::gpar(
                    fontsize = 7.5,
                    fontfamily = "Helvetica",
                    col = "grey15"
                )
            )
            grobs[[length(grobs) + 1L]] <- grid::textGrob(
                as.character(pw$p[[i]]),
                x = x_p,
                y = unit(yi, "npc"),
                just = c("left", "center"),
                gp = grid::gpar(
                    fontsize = 7.5,
                    fontfamily = "Helvetica",
                    col = "grey15"
                )
            )
        }
        inset <- do.call(grid::grobTree, grobs)
        box <- mc$pairwise_box
        fit_plot$plot <- fit_plot$plot +
            annotation_custom(
                grob = inset,
                xmin = box[["xmin"]],
                xmax = box[["xmax"]],
                ymin = box[["ymin"]],
                ymax = box[["ymax"]]
            )
    }

    fit_plot$plot <- fit_plot$plot +
        theme(
            panel.grid = element_blank(),
            legend.position = mc$legend_position
        )

    # ggsurvplot is not a ggplot; print() lays out plot + risk table.
    fmt <- tolower(tools::file_ext(io$output))
    out_dir <- dirname(io$output)
    if (!dir.exists(out_dir)) {
        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    }
    w <- mc$size$width
    h <- mc$size$height
    if (identical(fmt, "png")) {
        grDevices::png(
            io$output, width = w, height = h, units = "in",
            res = 300, bg = "white"
        )
    } else if (identical(fmt, "pdf")) {
        grDevices::pdf(io$output, width = w, height = h)
    } else if (identical(fmt, "svg")) {
        grDevices::svg(io$output, width = w, height = h, bg = "white")
    } else {
        stop(sprintf("Unsupported output format: %s", fmt), call. = FALSE)
    }
    tryCatch(
        print(fit_plot, newpage = FALSE),
        finally = grDevices::dev.off()
    )
    message("Saved: ", normalizePath(io$output, mustWork = TRUE))
    quit(save = "no")
}

# -----------------------------------------------------------------------------
# BASIC STYLE
# -----------------------------------------------------------------------------
pvalue <- round(surv_pvalue(fits, data = df)$pval, 3)
cox <- summary(coxph(f, data = df))
hr <- signif(cox$coef[, 2], digits = 2)
hr_lo <- signif(cox$conf.int[, "lower .95"], 3)
hr_up <- signif(cox$conf.int[, "upper .95"], 3)
ci_95 <- paste0(hr_lo, "-", hr_up)
# For >2 groups, cox$coef has multiple rows — keep first contrast text.
if (length(hr) > 1L) {
    hr <- hr[[1]]
    ci_95 <- paste0(hr_lo[[1]], "-", hr_up[[1]])
}
anno_label <- if (!is.null(config$pval_label) &&
    nzchar(as.character(config$pval_label))) {
    as.character(config$pval_label)
} else {
    sprintf("Log-rank p = %s\nHR = %s(%s)", pvalue, hr, ci_95)
}

group_table <- table(df[[group_col]])
legend_labels <- sprintf("%s (N = %s)", names(group_table), group_table)
palette_values <- unname(cols[seq_along(legend_labels)])

fit_plot <- ggsurvplot(
    fits,
    data = df,
    legend.title = "",
    legend = c(0.8, 0.9),
    title = config$labels$title,
    xlab = config$labels$x,
    ylab = config$labels$y,
    pval = FALSE,
    censor = TRUE,
    risk.table = FALSE,
    legend.labs = legend_labels,
    palette = palette_values,
    ggtheme = theme_prism() +
        theme(
            legend.title = element_text(
                face = "bold", family = "Times", colour = "black", size = 12
            ),
            legend.text = element_text(
                face = "bold", family = "Times", colour = "black", size = 12
            ),
            plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
            panel.border = element_rect(fill = "transparent", linewidth = 2),
            axis.line = element_blank(),
            panel.grid = element_blank()
        )
)

p <- fit_plot$plot +
    annotate(
        "text",
        x = 0,
        y = 0.2,
        hjust = 0,
        size = 4,
        label = anno_label
    ) +
    scale_x_continuous(expand = c(0.1, 0)) +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y
    ) +
    theme(plot.background = element_blank())

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output, width = config$size$width, height = config$size$height)
