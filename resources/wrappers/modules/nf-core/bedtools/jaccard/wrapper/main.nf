#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { BEDTOOLS_JACCARD } from '../main.nf'


params.a          = null
params.b          = null
params.outdir     = null

workflow {
    in_ch = Channel.value([[id: 'result', single_end: false], file(params.a, checkIfExists: true), file(params.b, checkIfExists: true)])

    BEDTOOLS_JACCARD(in_ch, [[], []])
}
