#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { PICARD_MERGESAMFILES } from '../main.nf'

params.bam1   = null
params.bam2   = null
params.outdir = null

workflow {
    bams_ch = Channel.value([
        [id: 'merged', single_end: false],
        [file(params.bam1, checkIfExists: true), file(params.bam2, checkIfExists: true)]
    ])

    PICARD_MERGESAMFILES(bams_ch)
}
