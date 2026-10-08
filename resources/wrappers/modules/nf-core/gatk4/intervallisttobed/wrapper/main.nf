#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { GATK4_INTERVALLISTTOBED } from '../main.nf'


params.intervals  = null
params.outdir     = null

workflow {
    intervals_ch = Channel.value([[id: 'test', single_end: false], file(params.intervals, checkIfExists: true)])

    GATK4_INTERVALLISTTOBED(intervals_ch)
}
