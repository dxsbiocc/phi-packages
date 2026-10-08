#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { CUSTOM_TABULARTOGSEACHIP } from '../main.nf'


params.tabular    = null
params.id_col     = 'gene_id'
params.symbol_col = 'gene_name'
params.outdir     = null

workflow {
    in_ch = Channel.value([[id: 'result'], file(params.tabular, checkIfExists: true)])

    CUSTOM_TABULARTOGSEACHIP(in_ch, Channel.value([params.id_col, params.symbol_col]))
}
