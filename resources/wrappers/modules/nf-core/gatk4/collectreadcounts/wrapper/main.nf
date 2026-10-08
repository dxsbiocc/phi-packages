#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { GATK4_COLLECTREADCOUNTS } from '../main.nf'


params.bam        = null
params.bai        = null
params.intervals  = null
params.outdir     = null

workflow {
    input_ch = Channel.value([[id: 'test', single_end: false], file(params.bam, checkIfExists: true), file(params.bai, checkIfExists: true), file(params.intervals, checkIfExists: true)])
    ref_ch   = Channel.value([[], []])

    GATK4_COLLECTREADCOUNTS(input_ch, ref_ch, ref_ch, ref_ch)
}
