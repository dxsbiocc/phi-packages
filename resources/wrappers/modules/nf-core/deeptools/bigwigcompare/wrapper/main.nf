#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { DEEPTOOLS_BIGWIGCOMPARE } from '../main.nf'

params.bigwig1 = null
params.bigwig2 = null
params.outdir  = null

workflow {
    bigwig_ch = Channel.value([[id: 'test'], file(params.bigwig1, checkIfExists: true), file(params.bigwig2, checkIfExists: true)])
    blacklist_ch = Channel.value([[id: 'no_blacklist'], []])

    DEEPTOOLS_BIGWIGCOMPARE(bigwig_ch, blacklist_ch)
}
