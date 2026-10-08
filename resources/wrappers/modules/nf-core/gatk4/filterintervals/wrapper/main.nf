#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { GATK4_FILTERINTERVALS } from '../main.nf'


params.intervals  = null
params.counts     = null
params.annotated  = null
params.outdir     = null

workflow {
    intervals_ch = Channel.value([[id: 'test'], file(params.intervals, checkIfExists: true)])
    counts_ch    = Channel.value([[:], [file(params.counts, checkIfExists: true)]])
    annotated_ch = Channel.value([[:], file(params.annotated, checkIfExists: true)])

    GATK4_FILTERINTERVALS(intervals_ch, counts_ch, annotated_ch)
}
