#!/usr/bin/env nextflow
// Thin agent-facing adapter over the DE_VOLCANO process at ../main.nf --
// takes a DE results table (e.g. from ../../deseq2/, ../../limma/, or
// ../../timeseries/ in this same family) and renders a volcano plot.
nextflow.enable.dsl = 2

include { DE_VOLCANO } from '../main.nf'

params.results       = null
params.gene_id_col   = 'gene_id'
params.log2fc_col    = 'log2FoldChange'
params.pvalue_col    = 'padj'
params.alpha         = '0.05'
params.lfc_threshold = '1'
params.label_top_n   = '10'
params.prefix        = 'de'
params.outdir        = null

workflow {
    results_ch = Channel.fromPath(params.results, checkIfExists: true)

    DE_VOLCANO(results_ch)
}
