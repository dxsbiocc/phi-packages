#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { CUSTOM_TABULARTOGSEACLS } from '../main.nf'


params.samples    = null
params.variable   = 'treatment'
params.outdir     = null

workflow {
    in_ch = Channel.value([[id: 'result'], file(params.samples, checkIfExists: true)])

    CUSTOM_TABULARTOGSEACLS(in_ch)
}
