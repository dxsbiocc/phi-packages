#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { UCSC_BIGWIGAVERAGEOVERBED } from '../main.nf'

params.bed    = null
params.bigwig = null
params.outdir = null

workflow {
    bed_ch    = Channel.value([[id: 'test'], [file(params.bed, checkIfExists: true)]])
    bigwig_ch = Channel.value(file(params.bigwig, checkIfExists: true))

    UCSC_BIGWIGAVERAGEOVERBED(bed_ch, bigwig_ch)
}
