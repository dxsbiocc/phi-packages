#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { CUSTOM_TABULARTOGSEAGCT } from '../main.nf'


params.tabular    = null
params.variable   = 'treatment'
params.reference  = 'mCherry'
params.target     = 'hND6'
params.outdir     = null

workflow {
    in_ch = Channel.value([[id: 'result', variable: params.variable, reference: params.reference, target: params.target, blocking: ''], file(params.tabular, checkIfExists: true)])

    CUSTOM_TABULARTOGSEAGCT(in_ch)
}
