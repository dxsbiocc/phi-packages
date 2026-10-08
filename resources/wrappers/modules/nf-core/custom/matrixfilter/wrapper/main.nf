#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { CUSTOM_MATRIXFILTER } from '../main.nf'


params.abundance  = null
params.samplesheet = null
params.outdir     = null

workflow {
    abundance_ch   = Channel.value([[id: 'result'], file(params.abundance, checkIfExists: true)])
    samplesheet_ch = Channel.value([[id: 'result'], file(params.samplesheet, checkIfExists: true)])

    CUSTOM_MATRIXFILTER(abundance_ch, samplesheet_ch)
}
