#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { BEDTOOLS_CLOSEST } from '../main.nf'

params.a      = null
params.b      = null
params.outdir = null

workflow {
    intervals_ch = Channel.value([[id: 'closest'], file(params.a, checkIfExists: true), file(params.b, checkIfExists: true)])
    fai_ch       = Channel.value([])

    BEDTOOLS_CLOSEST(intervals_ch, fai_ch)
}
