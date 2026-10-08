#!/usr/bin/env Rscript

# Template-ID: line-survival
#
# Purpose:
#   Draw a Kaplan-Meier survival curve by group, annotated with
#   log-rank p-value and Cox proportional-hazards HR.
#
# Inputs:
#   A table with one row per subject. Default example:
#     - OS_time: follow-up time
#     - OS: event indicator (1 = event, 0 = censored)
#     - group: strata
#
# Output:
#   A PDF, PNG, or SVG survival curve.
#
# Dependencies:
#   ggplot2, readr, survival, survminer, ggprism
#
# Example:
#   Rscript plot.R example.tsv output.pdf
#
# Agent adaptation:
#   For a new table, edit only CONFIG (column names and labels).
#   Edit DATA PREPARATION to change the Surv() formula or how HR
#   is formatted.
#   Edit PLOT only when the chart geometry or styling must change.
#
# Scientific assumptions:
#   Time and event columns define a right-censored survival object.
#   Groups are compared with a log-rank test and a Cox model.
#   Rows with missing time, event, or group are dropped.
#   The HR is for the second group relative to the first.

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
        time = "OS_time",
        event = "OS",
        group = "group"
    ),
    labels = list(
        title = "TCGA-BRCA",
        x = "Overall survival(Days)",
        y = "Survival Probability"
    )
)

load_packages(c("ggplot2", "readr", "survival", "survminer", "ggprism"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
df <- read_table_auto(io$input)
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
df[[group_col]] <- factor(df[[group_col]], levels = unique(df[[group_col]]))

f <- as.formula(sprintf(
    "Surv(%s, %s) ~ `%s`",
    time_col, event_col, group_col
))
fits <- survfit(f, data = df)
fits$call$formula <- f
pvalue <- round(surv_pvalue(fits, data = df)$pval, 3)
cox <- summary(coxph(f, data = df))
hr <- signif(cox$coef[, 2], digits = 2)
hr_lo <- signif(cox$conf.int[, "lower .95"], 3)
hr_up <- signif(cox$conf.int[, "upper .95"], 3)
ci_95 <- paste0(hr_lo, "-", hr_up)
anno_label <- sprintf("Log-rank p = %s\nHR = %s(%s)", pvalue, hr, ci_95)

group_table <- table(df[[group_col]])
legend_labels <- sprintf("%s (N = %s)", names(group_table), group_table)

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
palette_values <- c("#d46780", "#798234")

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
            axis.line = element_blank()
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
save_ggplot(p, io$output)
