#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { PICARD_SORTSAM } from '../main.nf'

params.bam        = null
params.sort_order = 'coordinate'
params.outdir     = null

workflow {
    bam_ch = Channel
        .fromPath(params.bam, checkIfExists: true)
        .map { bam -> [[id: "${bam.simpleName}_sorted", single_end: false], bam] }

    PICARD_SORTSAM(bam_ch, params.sort_order)
}
