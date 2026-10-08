#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { BEDTOOLS_SHUFFLE } from '../main.nf'


params.bed        = null
params.fai        = null
params.outdir     = null

workflow {
    bed_ch = Channel.value([[id: 'result'], file(params.bed, checkIfExists: true)])
    fai_ch = Channel.value([[id: 'genome'], file(params.fai, checkIfExists: true)])

    BEDTOOLS_SHUFFLE(bed_ch, fai_ch, [], [])
}
