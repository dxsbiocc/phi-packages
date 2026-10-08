#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { BEDTOOLS_SORT } from '../main.nf'

params.intervals = null
params.outdir    = null

workflow {
    intervals_ch = Channel.value([[id: 'sorted'], file(params.intervals, checkIfExists: true)])
    genome_ch    = Channel.value([])

    BEDTOOLS_SORT(intervals_ch, genome_ch)
}
