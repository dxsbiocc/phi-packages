#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { GATK4_GATHERBQSRREPORTS } from '../main.nf'


params.table      = null
params.outdir     = null

workflow {
    table_ch = Channel.value([[id: 'gathered', single_end: false], file(params.table, checkIfExists: true)])

    GATK4_GATHERBQSRREPORTS(table_ch)
}
