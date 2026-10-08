#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { PICARD_COLLECTINSERTSIZEMETRICS } from '../main.nf'

params.bam    = null
params.outdir = null

workflow {
    bam_ch = Channel
        .fromPath(params.bam, checkIfExists: true)
        .map { bam -> [[id: bam.simpleName, single_end: false], bam] }

    PICARD_COLLECTINSERTSIZEMETRICS(bam_ch)
}
