#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { DEEPTOOLS_MULTIBIGWIGSUMMARY } from '../main.nf'

params.bigwig1 = null
params.bigwig2 = null
params.outdir  = null

workflow {
    input_ch = Channel.value([
        [id: 'test'],
        [file(params.bigwig1, checkIfExists: true), file(params.bigwig2, checkIfExists: true)],
        ['bigwig1', 'bigwig2']
    ])
    blacklist_ch = Channel.value([[id: 'no_blacklist'], []])

    DEEPTOOLS_MULTIBIGWIGSUMMARY(input_ch, blacklist_ch)
}
