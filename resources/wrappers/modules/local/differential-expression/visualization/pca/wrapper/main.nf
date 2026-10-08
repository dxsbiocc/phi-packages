#!/usr/bin/env nextflow
// Thin agent-facing adapter over the DE_PCA process at ../main.nf -- takes
// a variance-stabilized/rlog expression matrix (e.g. from ../../deseq2/'s
// vst_counts output in this same family) plus a sample sheet, and renders
// a PCA scatter colored (and optionally shaped) by sample metadata.
nextflow.enable.dsl = 2

include { DE_PCA } from '../main.nf'

params.matrix        = null
params.samplesheet   = null
params.gene_id_col   = 'gene_id'
params.sample_id_col = 'sample'
params.color_by      = 'condition'
params.shape_by      = 'null'
params.ntop          = '500'
params.prefix        = 'de'
params.outdir        = null

workflow {
    matrix_ch      = Channel.fromPath(params.matrix, checkIfExists: true)
    samplesheet_ch = Channel.fromPath(params.samplesheet, checkIfExists: true)

    DE_PCA(matrix_ch, samplesheet_ch)
}
