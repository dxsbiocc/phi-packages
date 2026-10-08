#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { DEEPTOOLS_BAMCOVERAGE } from '../main.nf'

params.bam    = null
params.bai    = null
params.outdir = null

workflow {
    bam_ch = Channel.value([[id: 'test', single_end: false], file(params.bam, checkIfExists: true), file(params.bai, checkIfExists: true)])
    empty_ch = Channel.value([])
    blacklist_ch = Channel.value([[id: 'no_blacklist'], []])

    DEEPTOOLS_BAMCOVERAGE(bam_ch, empty_ch, empty_ch, blacklist_ch)
}
