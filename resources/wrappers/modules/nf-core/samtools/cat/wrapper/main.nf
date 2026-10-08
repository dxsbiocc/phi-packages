#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { SAMTOOLS_CAT } from '../main.nf'

params.bam1   = null
params.bam2   = null
params.outdir = null

workflow {
    bams_ch = Channel.value([[id: 'test'], [file(params.bam1, checkIfExists: true), file(params.bam2, checkIfExists: true)]])

    SAMTOOLS_CAT(bams_ch)
}
