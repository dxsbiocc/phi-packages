#!/usr/bin/env Rscript

# Classification logic adapted from bioSkills's differential-expression/
# de-visualization volcano_plot.R example (MIT license, GPTomics).
# Styling (diverging direction colors, point-size-by-magnitude, restrained
# theme) adapted from resources/skills/omics-visualization's scatter/volcano
# template -- values and layout only, not sourced code, so this script
# stays self-contained and runnable in its own container. See
# resources/palettes/palettes.yaml
# (Diverging.RdBu) for where the direction colors come from.

################################################
################################################
## Functions                                  ##
################################################
################################################

#' Parse out options from a string without recourse to optparse
#'
#' @param x Long-form argument list like --opt1 val1 --opt2 val2
#' @return named list of options and values similar to optparse
parse_args <- function(x) {
    args_list <- unlist(strsplit(x, ' ?--')[[1]])[-1]
    args_vals <- lapply(args_list, function(x) scan(text = x, what = 'character', quiet = TRUE))
    args_vals <- lapply(args_vals, function(z) { length(z) <- 2; z })
    parsed_args <- structure(lapply(args_vals, function(x) x[2]), names = lapply(args_vals, function(x) x[1]))
    parsed_args[!is.na(parsed_args)]
}

#' Flexibly read CSV or TSV files
read_delim_flexible <- function(file, header = TRUE, row.names = NULL, check.names = TRUE) {
    ext <- tolower(tail(strsplit(basename(file), split = "\\\\.")[[1]], 1))
    if (ext == "tsv" || ext == "txt") {
        separator <- "\\t"
    } else if (ext == "csv") {
        separator <- ","
    } else {
        stop(paste("Unknown separator for", ext))
    }
    read.delim(file, sep = separator, header = header, row.names = row.names, check.names = check.names)
}

################################################
################################################
## PARSE PARAMETERS FROM NEXTFLOW             ##
################################################
################################################

opt <- list(
    output_prefix = ifelse('$task.ext.prefix' == 'null', 'de', '$task.ext.prefix'),
    results_file  = '$results',
    gene_id_col   = 'gene_id',
    log2fc_col    = 'log2FoldChange',
    pvalue_col    = 'padj',
    alpha         = 0.05,
    lfc_threshold = 1,
    label_top_n   = 10
)
opt_types <- lapply(opt, class)

args_opt <- parse_args('$task.ext.args')
for (ao in names(args_opt)) {
    if (!ao %in% names(opt)) stop(paste("Invalid option:", ao))
    opt[[ao]] <- as(args_opt[[ao]], opt_types[[ao]])
}

################################################
################################################
## LOAD LIBRARIES                             ##
################################################
################################################

library(ggplot2)

################################################
################################################
## READ AND CLASSIFY RESULTS                  ##
################################################
################################################

res <- read_delim_flexible(opt\$results_file, header = TRUE, check.names = FALSE)

for (needed in c(opt\$gene_id_col, opt\$log2fc_col, opt\$pvalue_col)) {
    if (!needed %in% colnames(res)) {
        stop(paste0("Column '", needed, "' not found. Available columns: ", paste(colnames(res), collapse = ", ")))
    }
}

res <- res[stats::complete.cases(res[, c(opt\$log2fc_col, opt\$pvalue_col)]), ]
res\$log2fc      <- as.numeric(res[[opt\$log2fc_col]])
res\$qvalue       <- pmax(as.numeric(res[[opt\$pvalue_col]]), 1e-30) # floored so -log10 stays finite
res\$neg_log10_q  <- -log10(res\$qvalue)
res\$gene         <- res[[opt\$gene_id_col]]

# padj < alpha AND |log2FC| > lfc_threshold: standard significance call --
# see bioSkills's deseq2-basics/de-results skills for the tradeoffs of this
# cutoff versus TREAT/lfcThreshold-based FDR control.
res\$direction <- factor(
    ifelse(
        res\$qvalue < opt\$alpha & abs(res\$log2fc) > opt\$lfc_threshold,
        ifelse(res\$log2fc > 0, "Up", "Down"),
        "None"
    ),
    levels = c("Down", "None", "Up")
)

top_labels <- do.call(rbind, lapply(
    split(res, res\$direction),
    function(sub) {
        if (unique(sub\$direction) == "None" || nrow(sub) == 0) return(sub[0, ])
        sub[head(order(sub\$qvalue), opt\$label_top_n), ]
    }
))

################################################
################################################
## PLOT                                       ##
################################################
################################################

# Direction colors from Phi's Diverging.RdBu catalog entry
# (resources/palettes/palettes.yaml) -- the conventional red-up/blue-down reading
# for signed fold changes.
direction_colors <- c(Down = "#2166ac", None = "grey70", Up = "#b2182b")

PlotFile <- paste0(opt\$output_prefix, ".volcano.pdf")

p <- ggplot(res, aes(x = log2fc, y = neg_log10_q, colour = direction, label = gene)) +
    geom_point(aes(size = abs(log2fc)), alpha = 0.7, stroke = 0) +
    geom_text(data = top_labels, colour = "black", size = 3, vjust = -0.6, check_overlap = TRUE) +
    geom_vline(xintercept = c(-opt\$lfc_threshold, opt\$lfc_threshold), linetype = "dashed", colour = "grey40") +
    geom_hline(yintercept = -log10(opt\$alpha), linetype = "dashed", colour = "grey40") +
    scale_colour_manual(values = direction_colors, drop = FALSE) +
    scale_size_continuous(range = c(0.6, 4), guide = "none") +
    labs(
        title = "Volcano plot",
        subtitle = paste0("padj < ", opt\$alpha, ", |log2FC| > ", opt\$lfc_threshold),
        x = paste0(opt\$log2fc_col, " (log2 fold change)"),
        y = paste0("-log10(", opt\$pvalue_col, ")"),
        colour = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
        legend.position = "top",
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "grey40", fill = NA, linewidth = 0.5)
    )

ggsave(PlotFile, p, width = 7, height = 7, dpi = 300)

write.table(
    top_labels[, c("gene", "log2fc", "qvalue", "direction")],
    file = paste0(opt\$output_prefix, ".volcano.labels.tsv"),
    row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE
)

################################################
################################################
## VERSIONS FILE                              ##
################################################
################################################

r.version <- paste(R.version[['major']], R.version[['minor']], sep = ".")

writeLines(
    c(
        '"${task.process}":',
        paste('    r-base:', r.version),
        paste('    r-ggplot2:', as.character(packageVersion('ggplot2')))
    ),
    'versions.yml')
