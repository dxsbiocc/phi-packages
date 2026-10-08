#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { DEEPTOOLS_MULTIBAMSUMMARY } from '../main.nf'

params.bam1   = null
params.bai1   = null
params.bam2   = null
params.bai2   = null
params.outdir = null

workflow {
    input_ch = Channel.value([
        [id: 'test', single_end: false],
        [file(params.bam1, checkIfExists: true), file(params.bam2, checkIfExists: true)],
        [file(params.bai1, checkIfExists: true), file(params.bai2, checkIfExists: true)],
        ['bam1', 'bam2']
    ])
    blacklist_ch = Channel.value([[id: 'no_blacklist'], []])

    DEEPTOOLS_MULTIBAMSUMMARY(input_ch, blacklist_ch)
}
