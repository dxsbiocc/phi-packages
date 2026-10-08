#!/usr/bin/env Rscript

# Template-ID: boxplot-raincloud-differential
#
# Purpose:
#   Draw a flipped two-group raincloud: opposite half-violins, jitter,
#   quantile connectors colored by quantile-regression significance,
#   and KS / Wilcoxon annotations.
#
# Inputs:
#   A two-column table with exactly two groups. Default example:
#     - gender: group (two levels)
#     - TP53: numeric value
#
# Output:
#   A PDF, PNG, or SVG differential raincloud plot.
#
# Dependencies:
#   ggplot2, readr, gghalves, quantreg, ggprism
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (column names and labels).
#   Edit DATA PREPARATION to change group order, quantile grid, or
#   annotation positions.
#   Edit PLOT only when the chart geometry or styling must change.
#
# Scientific assumptions:
#   The x column has exactly two groups.
#   Quantile lines are from a linear quantile regression at
#   tau = 0.1, 0.2, ..., 0.9; color marks p < 0.05 on the group
#   coefficient (bootstrap SE, R = 1000).
#   KS and Wilcoxon tests compare the two full distributions.
#   Flipping axes is a display choice.

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
        x = "gender",
        y = "TP53"
    ),
    labels = list(
        title = "",
        x = "Group",
        y = "Gene expression level (TP53)"
    )
)

load_packages(c("ggplot2", "readr", "gghalves", "quantreg", "ggprism"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
require_columns(df, config$columns)
x_col <- config$columns$x
y_col <- config$columns$y
df[[y_col]] <- as.numeric(df[[y_col]])
x_names <- sort(unique(as.character(df[[x_col]])))
if (length(x_names) != 2L) {
    stop("This plot requires exactly two groups in the x column.", call. = FALSE)
}
df[[x_col]] <- factor(df[[x_col]], levels = x_names)

hjust <- 0.05
df$x <- ifelse(df[[x_col]] == x_names[[1]], 1, 2)
df$vx <- ifelse(
    df[[x_col]] == x_names[[1]],
    1 - hjust,
    2 + hjust
)
left <- df[df[[x_col]] == x_names[[1]], ]
right <- df[df[[x_col]] == x_names[[2]], ]

f <- as.formula(paste(y_col, "~", x_col))
suppressWarnings(
    fit <- rq(f, tau = seq(0.1, 0.9, 0.1), data = df)
)
pvalues <- unlist(lapply(
    summary(fit, se = "boot", R = 1000),
    function(x) x$coefficients[2, "Pr(>|t|)"]
))

probs <- seq(0.1, 0.9, by = 0.1)
qual_rows <- list()
for (g in x_names) {
    vals <- df[[y_col]][df[[x_col]] == g]
    qual_rows[[g]] <- data.frame(
        x = ifelse(g == x_names[[1]], 1 + hjust * 2, 2 - hjust * 2),
        q = as.numeric(quantile(vals, probs, na.rm = TRUE)),
        prob = probs,
        signif = ifelse(pvalues < 0.05, "P < 0.05", "NS")
    )
    qual_rows[[g]][[x_col]] <- g
}
qual_data <- do.call(rbind, qual_rows)

format_p_value <- function(p) {
    if (p < 0.01) {
        formatC(p, format = "e", digits = 2)
    } else {
        round(p, 4)
    }
}
grouped <- split(df[[y_col]], df[[x_col]])
suppressWarnings(
    ks_pvalue <- ks.test(grouped[[1]], grouped[[2]])[["p.value"]]
)
suppressWarnings(
    wt_pvalue <- wilcox.test(grouped[[1]], grouped[[2]])[["p.value"]]
)
anno_labels <- c(
    paste0("Kolmogorov-Smirnov Test\nP = ", format_p_value(ks_pvalue)),
    paste0("Wilcoxon Rank Sum Test\nP = ", format_p_value(wt_pvalue))
)

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
fill_values <- c("#a0c443", "#5da4dc")

p <- ggplot(df, aes(y = .data[[y_col]])) +
    geom_half_violin(
        data = left,
        aes(x = vx, fill = .data[[x_col]]),
        side = "l",
        colour = "grey50"
    ) +
    geom_half_violin(
        data = right,
        aes(x = vx, fill = .data[[x_col]]),
        side = "r",
        colour = "grey50"
    ) +
    geom_jitter(aes(x = x), width = 0.03) +
    geom_line(
        data = qual_data,
        aes(x = x, y = q, group = prob, colour = signif)
    ) +
    geom_point(
        data = qual_data,
        aes(x = x, y = q, colour = signif),
        size = 4
    ) +
    annotate(
        "text",
        x = c(1.3, 1.6),
        y = 4.3,
        label = anno_labels
    ) +
    scale_fill_manual(values = fill_values) +
    scale_colour_manual(
        values = c("NS" = "grey80", "P < 0.05" = "black")
    ) +
    scale_x_continuous(breaks = c(1, 2), labels = x_names) +
    coord_flip() +
    labs(
        title = config$labels$title,
        x = config$labels$x,
        y = config$labels$y
    ) +
    theme_prism() +
    theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.background = element_blank()
    )

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------
save_ggplot(p, io$output)
