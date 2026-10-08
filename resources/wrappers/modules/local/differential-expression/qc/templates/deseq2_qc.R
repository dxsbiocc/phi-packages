#!/usr/bin/env Rscript

# Adapted from nf-core/rnaseq's bin/deseq2_qc.r (Harshil Patel, Gavin Kelly;
# MIT license) -- rewritten from optparse/CLI-flag parsing to Nextflow's
# template variable interpolation, matching the other scripts in this
# differential-expression family (deseq2_differential.R, limma_de.R,
# dream.R). The original script's MultiQC section-header composition step
# is dropped here -- this wrapper's outputs never surfaced those files.

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

##' PCA pre-processor
##'
##' Generate all the necessary information to plot PCA from a DESeq2 object
##' in which an assay containing a variance-stabilised matrix of counts is
##' stored. Copied from DESeq2::plotPCA, but with additional ability to
##' say which assay to run the PCA on.
##'
##' @param object The DESeq2DataSet object.
##' @param ntop number of top genes to use for principal components, selected by highest row variance.
##' @param assay the name or index of the assay that stores the variance-stabilised data.
##' @return A data.frame containing the projected data alongside the grouping columns.
##' A 'percentVar' attribute is set which includes the percentage of variation each PC explains.
##' @author Gavin Kelly
plotPCA_vst <- function(object, ntop = 500, assay = length(assays(object))) {
    rv         <- rowVars(assay(object, assay))
    select     <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
    pca        <- prcomp(t(assay(object, assay)[select, ]), center = TRUE, scale = FALSE)
    percentVar <- pca\$sdev^2 / sum(pca\$sdev^2)
    df         <- cbind(as.data.frame(colData(object)), pca\$x)
    # Order points so extreme samples are more likely to get a label.
    ord        <- order(abs(rank(df\$PC1) - median(df\$PC1)), abs(rank(df\$PC2) - median(df\$PC2)))
    df         <- df[ord, ]
    attr(df, "percentVar") <- data.frame(PC = seq(along = percentVar), percentVar = 100 * percentVar)
    return(df)
}

################################################
################################################
## PARSE PARAMETERS FROM NEXTFLOW             ##
################################################
################################################

opt <- list(
    output_prefix = ifelse('$task.ext.prefix' == 'null', 'deseq2', '$task.ext.prefix'),
    count_file    = '$counts',
    count_col     = 3,
    id_col        = 1,
    sample_suffix = '',
    vst           = FALSE
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

library(DESeq2)
library(ggplot2)
library(pheatmap)

# Palette hexes copied from Phi's resources/palettes/palettes.yaml
# (Quantitative.BluGrn, the shared heatmap
# default, and a legible dark hue from Qualitative.Safe for points/text on
# white) -- values only, not code, so this script stays self-contained and
# runnable in its own container.
PALETTE_HEATMAP_SEQUENTIAL <- c("#c4e6c3", "#96d2a4", "#6dbc90", "#4da284", "#36877a", "#266b6e", "#1d4f60")
PALETTE_POINT <- "#332288"

################################################
################################################
## READ IN COUNTS FILE                        ##
################################################
################################################

count.table            <- read.delim(file = opt\$count_file, header = TRUE, row.names = NULL, check.names = FALSE)
rownames(count.table)  <- count.table[, opt\$id_col]
count.table            <- count.table[, opt\$count_col:ncol(count.table), drop = FALSE]
colnames(count.table)  <- gsub(opt\$sample_suffix, "", colnames(count.table))
colnames(count.table)  <- gsub(pattern = '\\\\.\$', replacement = '', colnames(count.table))

################################################
################################################
## RUN DESEQ2                                 ##
################################################
################################################

samples.vec     <- colnames(count.table)
name_components <- strsplit(samples.vec, "_")
n_components    <- length(name_components[[1]])
decompose       <- n_components != 1 && all(sapply(name_components, length) == n_components)
coldata         <- data.frame(sample = samples.vec, row.names = samples.vec)
if (decompose) {
    groupings        <- as.data.frame(lapply(1:n_components, function(i) sapply(name_components, "[[", i)))
    n_distinct       <- sapply(groupings, function(grp) length(unique(grp)))
    groupings        <- groupings[n_distinct != 1 & n_distinct != length(samples.vec)]
    if (ncol(groupings) != 0) {
        names(groupings) <- paste0("Group", 1:ncol(groupings))
        coldata <- cbind(coldata, groupings)
    } else {
        decompose <- FALSE
    }
}

DDSFile <- paste(opt\$output_prefix, ".dds.RData", sep = "")

counts <- count.table[, samples.vec, drop = FALSE]
# `design=~1` creates intercept-only model, equivalent to setting `blind=TRUE` for transformation.
dds    <- DESeqDataSetFromMatrix(countData = round(counts), colData = coldata, design = ~1)
dds    <- estimateSizeFactors(dds)
if (min(dim(count.table)) <= 1) { # No point if only one sample, or one gene
    save(dds, file = DDSFile)
    saveRDS(dds, file = sub("\\\\.dds\\\\.RData\$", ".rds", DDSFile))
    warning("Not enough samples or genes in counts file for PCA.", call. = FALSE)
    quit(save = "no", status = 0, runLast = FALSE)
}
if (!opt\$vst) {
    vst_name <- "rlog"
    rld      <- rlog(dds)
} else {
    vst_name <- "vst"
    rld      <- varianceStabilizingTransformation(dds)
}

assay(dds, vst_name) <- assay(rld)
save(dds, file = DDSFile)
saveRDS(dds, file = sub("\\\\.dds\\\\.RData\$", ".rds", DDSFile))

################################################
################################################
## PLOT QC                                    ##
################################################
################################################

PlotFile <- paste(opt\$output_prefix, ".plots.pdf", sep = "")

pdf(file = PlotFile, onefile = TRUE, width = 7, height = 7)
## PCA
ntop <- c(500, Inf)
for (n_top_var in ntop) {
    pca.data      <- plotPCA_vst(dds, assay = vst_name, ntop = n_top_var)
    percentVar    <- round(attr(pca.data, "percentVar")\$percentVar)
    plot_subtitle <- ifelse(n_top_var == Inf, "All genes", paste("Top", n_top_var, "genes"))
    pl <- ggplot(pca.data, aes(PC1, PC2, label = paste0(" ", sample, " "))) +
        geom_point(color = PALETTE_POINT, size = 2) +
        geom_text(check_overlap = TRUE, vjust = 0.5, hjust = "inward", color = PALETTE_POINT, size = 3.2) +
        xlab(paste0("PC1: ", percentVar[1], "% variance")) +
        ylab(paste0("PC2: ", percentVar[2], "% variance")) +
        labs(title = paste0("First PCs on ", vst_name, "-transformed data"), subtitle = plot_subtitle) +
        theme_minimal(base_size = 11) +
        theme(legend.position = "top",
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "grey40", fill = NA, linewidth = 0.5))
    print(pl)

    if (decompose) {
        pc_names    <- paste0("PC", attr(pca.data, "percentVar")\$PC)
        long_pc     <- reshape(pca.data, varying = pc_names, direction = "long", sep = "", timevar = "component", idvar = "pcrow")
        long_pc     <- subset(long_pc, component <= 5)
        long_pc_grp <- reshape(long_pc, varying = names(groupings), direction = "long", sep = "", timevar = "grouper")
        long_pc_grp <- subset(long_pc_grp, grouper <= 5)
        long_pc_grp\$component <- paste("PC", long_pc_grp\$component)
        long_pc_grp\$grouper   <- paste0(long_pc_grp\$grouper, c("st", "nd", "rd", "th", "th")[long_pc_grp\$grouper], " prefix")
        pl <- ggplot(long_pc_grp, aes(x = Group, y = PC)) +
            geom_point() +
            stat_summary(fun = mean, geom = "line", aes(group = 1)) +
            labs(x = NULL, y = NULL, subtitle = plot_subtitle, title = "PCs split by sample-name prefixes") +
            facet_grid(component ~ grouper, scales = "free_x") +
            scale_x_discrete(guide = guide_axis(n.dodge = 3))
        print(pl)
    }
} # at end of loop, we'll be using the user-defined ntop if any, else all genes

## WRITE PC1 vs PC2 VALUES TO FILE
pca.vals           <- pca.data[, c("PC1", "PC2")]
colnames(pca.vals) <- paste0(colnames(pca.vals), ": ", percentVar[1:2], '% variance')
pca.vals           <- cbind(sample = rownames(pca.vals), pca.vals)
write.table(pca.vals, file = paste(opt\$output_prefix, ".pca.vals.txt", sep = ""),
            row.names = FALSE, col.names = TRUE, sep = "\t", quote = TRUE)

## SAMPLE CORRELATION HEATMAP
sampleDists      <- dist(t(assay(dds, vst_name)))
sampleDistMatrix <- as.matrix(sampleDists)
colors           <- colorRampPalette(rev(PALETTE_HEATMAP_SEQUENTIAL))(255)
pheatmap(
    sampleDistMatrix,
    clustering_distance_rows = sampleDists,
    clustering_distance_cols = sampleDists,
    col  = colors,
    main = paste("Euclidean distance between", vst_name, "of samples")
)

## WRITE SAMPLE DISTANCES TO FILE
write.table(cbind(sample = rownames(sampleDistMatrix), sampleDistMatrix), file = paste(opt\$output_prefix, ".sample.dists.txt", sep = ""),
            row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
dev.off()

################################################
## OPTIONAL MULTIQC FRAGMENTS                 ##
################################################

label_lower <- tolower('$task.ext.args2')
if (nzchar(label_lower) && label_lower != 'null' &&
    file.exists('$pca_header_multiqc') && file.exists('$clustering_header_multiqc')) {
    pca_values <- paste(opt\$output_prefix, ".pca.vals.txt", sep = "")
    distance_values <- paste(opt\$output_prefix, ".sample.dists.txt", sep = "")
    writeLines(
        c(readLines('$pca_header_multiqc'), readLines(pca_values)),
        paste0(label_lower, ".pca.vals_mqc.tsv")
    )
    writeLines(
        c(readLines('$clustering_header_multiqc'), readLines(distance_values)),
        paste0(label_lower, ".sample.dists_mqc.tsv")
    )
}

################################################
################################################
## SAVE SIZE FACTORS                          ##
################################################
################################################

SizeFactorsDir <- "size_factors/"
if (file.exists(SizeFactorsDir) == FALSE) {
    dir.create(SizeFactorsDir, recursive = TRUE)
}

NormFactorsFile <- paste(SizeFactorsDir, opt\$output_prefix, ".size_factors.RData", sep = "")

normFactors <- sizeFactors(dds)
save(normFactors, file = NormFactorsFile)

for (name in names(sizeFactors(dds))) {
    sizeFactorFile <- paste(SizeFactorsDir, name, ".txt", sep = "")
    write(as.numeric(sizeFactors(dds)[name]), file = sizeFactorFile)
}

################################################
################################################
## R SESSION INFO                             ##
################################################
################################################

sink("R_sessionInfo.log")
print(sessionInfo())
sink()

################################################
################################################
## VERSIONS FILE                              ##
################################################
################################################

r.version      <- paste(R.version[['major']], R.version[['minor']], sep = ".")
deseq2.version <- as.character(packageVersion('DESeq2'))

writeLines(
    c(
        '"${task.process}":',
        paste('    r-base:', r.version),
        paste('    bioconductor-deseq2:', deseq2.version)
    ),
    'versions.yml')
