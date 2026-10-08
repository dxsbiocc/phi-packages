#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { RTGTOOLS_SVDECOMPOSE } from '../main.nf'


params.sv         = null
params.sv_tbi     = null
params.outdir     = null

workflow {
    in_ch = Channel.value([[id: 'sample'], file(params.sv, checkIfExists: true), file(params.sv_tbi, checkIfExists: true)])

    RTGTOOLS_SVDECOMPOSE(in_ch)
}
