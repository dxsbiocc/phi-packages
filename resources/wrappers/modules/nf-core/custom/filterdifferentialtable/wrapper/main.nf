#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { CUSTOM_FILTERDIFFERENTIALTABLE } from '../main.nf'


params.table      = null
params.logfc_column = 'log2FoldChange'
params.fc_threshold = 2
params.fc_cardinality = '>='
params.stat_column = 'padj'
params.stat_threshold = '0.05'
params.stat_cardinality = '<='
params.outdir     = null

workflow {
    in_ch = Channel.value([[id: 'result'], file(params.table, checkIfExists: true)])
    fc_ch   = Channel.value([params.logfc_column, params.fc_threshold, params.fc_cardinality])
    stat_ch = Channel.value([params.stat_column, params.stat_threshold, params.stat_cardinality])

    CUSTOM_FILTERDIFFERENTIALTABLE(in_ch, fc_ch, stat_ch)
}
