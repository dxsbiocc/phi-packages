#!/usr/bin/env nextflow
// Thin agent-facing adapter over the DE_HEATMAP process at ../main.nf --
// takes an expression matrix (e.g. from ../../deseq2/'s vst_counts output
// in this same family) plus a sample sheet, and renders a clustered
// heatmap of the top-variable genes, row z-scored.
nextflow.enable.dsl = 2

include { DE_HEATMAP } from '../main.nf'

params.matrix        = null
params.samplesheet   = null
params.gene_id_col   = 'gene_id'
params.sample_id_col = 'sample'
params.annotate_by   = 'condition'
params.top_n         = '50'
params.scale_rows    = true
params.prefix        = 'de'
params.outdir        = null

workflow {
    matrix_ch      = Channel.fromPath(params.matrix, checkIfExists: true)
    samplesheet_ch = Channel.fromPath(params.samplesheet, checkIfExists: true)

    DE_HEATMAP(matrix_ch, samplesheet_ch)
}
