#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { DEEPTOOLS_PLOTPROFILE } from '../main.nf'

params.matrix = null
params.outdir = null

workflow {
    matrix_ch = Channel.value([[id: 'test', single_end: false], file(params.matrix, checkIfExists: true)])

    DEEPTOOLS_PLOTPROFILE(matrix_ch)
}
