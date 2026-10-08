#!/usr/bin/env Rscript

# PCA computation adapted from bioSkills's differential-expression/
# de-visualization pca_plot.R example (MIT license, GPTomics) and this
# family's own ../../qc/templates/deseq2_qc.R (plotPCA_vst helper).
# Styling (colorblind-safe categorical palette, restrained theme) adapted
# from resources/skills/omics-visualization's scatter/group template --
# values and layout only, not sourced code, so this script stays
# self-contained and runnable in its own container. See
# resources/palettes/palettes.yaml
# (Qualitative.Safe) for where the point colors come from.

################################################
################################################
## Functions                                  ##
################################################
################################################

#' Parse out options from a string without recourse to optparse
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

#' Turn "null" or empty strings into actual NULL
nullify <- function(x) {
    if (is.character(x) && (tolower(x) == "null" || x == "")) NULL else x
}

################################################
################################################
## PARSE PARAMETERS FROM NEXTFLOW             ##
################################################
################################################

opt <- list(
    output_prefix = ifelse('$task.ext.prefix' == 'null', 'de', '$task.ext.prefix'),
    matrix_file    = '$matrix',
    sample_file    = '$samplesheet',
    gene_id_col    = 'gene_id',
    sample_id_col  = 'sample',
    color_by       = 'condition',
    shape_by       = "null",
    ntop           = 500
)
opt_types <- lapply(opt, class)

args_opt <- parse_args('$task.ext.args')
for (ao in names(args_opt)) {
    if (!ao %in% names(opt)) stop(paste("Invalid option:", ao))
    opt[[ao]] <- as(args_opt[[ao]], opt_types[[ao]])
}
opt\$shape_by <- nullify(opt\$shape_by)

################################################
################################################
## LOAD LIBRARIES                             ##
################################################
################################################

library(ggplot2)

################################################
################################################
## READ INPUTS                                ##
################################################
################################################

mat.table              <- read_delim_flexible(opt\$matrix_file, header = TRUE, row.names = NULL, check.names = FALSE)
rownames(mat.table)     <- mat.table[[opt\$gene_id_col]]
mat.table[[opt\$gene_id_col]] <- NULL
mat                     <- as.matrix(mat.table)

sample.sheet <- read_delim_flexible(opt\$sample_file, header = TRUE, check.names = FALSE)
if (!opt\$sample_id_col %in% colnames(sample.sheet)) {
    stop(paste0("Specified sample ID column '", opt\$sample_id_col, "' is not in the sample sheet"))
}
rownames(sample.sheet) <- sample.sheet[[opt\$sample_id_col]]

missing_samples <- setdiff(colnames(mat), rownames(sample.sheet))
if (length(missing_samples) > 0) {
    stop(paste("Samples in matrix missing from sample sheet:", paste(missing_samples, collapse = ",")))
}
sample.sheet <- sample.sheet[colnames(mat), , drop = FALSE]

for (needed in c(opt\$color_by, opt\$shape_by)) {
    if (!is.null(needed) && !needed %in% colnames(sample.sheet)) {
        stop(paste0("Column '", needed, "' not found in sample sheet. Available columns: ", paste(colnames(sample.sheet), collapse = ", ")))
    }
}

################################################
################################################
## RUN PCA                                    ##
################################################
################################################

##' PCA pre-processor -- selects the top-variance rows and runs prcomp(),
##' matching DESeq2::plotPCA's own approach (see ../../qc/templates/deseq2_qc.R).
run_pca <- function(mat, ntop) {
    rv     <- apply(mat, 1, stats::var) # base R only -- no matrixStats dependency
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
    pca    <- prcomp(t(mat[select, , drop = FALSE]), center = TRUE, scale. = FALSE)
    percentVar <- pca\$sdev^2 / sum(pca\$sdev^2)
    list(scores = as.data.frame(pca\$x), percentVar = percentVar)
}

pca_result       <- run_pca(mat, opt\$ntop)
pca.data         <- cbind(pca_result\$scores, sample.sheet) # rownames(sample.sheet) already match rownames(scores)
pca.data\$.sample <- rownames(pca.data)                     # own column, distinct from any samplesheet column
percentVar       <- round(100 * pca_result\$percentVar)

################################################
################################################
## PLOT                                       ##
################################################
################################################

# Colorblind-safe categorical hues from Phi's
# Qualitative.Safe catalog entry (resources/palettes/palettes.yaml).
QUALITATIVE_SAFE <- c("#88ccee", "#cc6677", "#ddcc77", "#117733", "#332288", "#aa4499",
                       "#44aa99", "#999933", "#882255", "#661100", "#6699cc", "#888888")

PlotFile <- paste0(opt\$output_prefix, ".pca.pdf")

aes_args <- list(x = quote(PC1), y = quote(PC2), colour = as.name(opt\$color_by))
if (!is.null(opt\$shape_by)) aes_args\$shape <- as.name(opt\$shape_by)

p <- ggplot(pca.data, do.call(aes, aes_args)) +
    geom_point(size = 3, alpha = 0.9) +
    scale_colour_manual(values = QUALITATIVE_SAFE) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    labs(title = "PCA of samples", subtitle = paste("Top", opt\$ntop, "variable genes")) +
    theme_minimal(base_size = 11) +
    theme(
        legend.position = "right",
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "grey40", fill = NA, linewidth = 0.5)
    )

ggsave(PlotFile, p, width = 7, height = 6, dpi = 300)

out.vals <- pca.data[, c(".sample", "PC1", "PC2")]
colnames(out.vals)[1] <- "sample"
write.table(
    out.vals,
    file = paste0(opt\$output_prefix, ".pca.vals.tsv"),
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
