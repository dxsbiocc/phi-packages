#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { RTGTOOLS_PEDFILTER } from '../main.nf'


params.ped        = null
params.outdir     = null

workflow {
    in_ch = Channel.value([[id: 'sample', single_end: false], file(params.ped, checkIfExists: true)])

    RTGTOOLS_PEDFILTER(in_ch)
}
