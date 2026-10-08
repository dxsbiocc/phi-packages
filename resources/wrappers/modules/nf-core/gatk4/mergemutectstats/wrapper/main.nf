#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { GATK4_MERGEMUTECTSTATS } from '../main.nf'


params.stats      = null
params.outdir     = null

workflow {
    stats_ch = Channel.value([[id: 'merged', single_end: false], file(params.stats, checkIfExists: true)])

    GATK4_MERGEMUTECTSTATS(stats_ch)
}
