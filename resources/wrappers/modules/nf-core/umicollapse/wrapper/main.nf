#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { UMICOLLAPSE } from '../main.nf'

params.bam        = null
params.bai        = null
params.single_end = false
params.mode       = 'bam'
params.outdir     = null

workflow {
    bam_bai_ch = Channel
        .fromPath(params.bam, checkIfExists: true)
        .map { bam -> [[id: bam.simpleName, single_end: params.single_end], bam] }
        .combine(Channel.fromPath(params.bai, checkIfExists: true))
        .map { meta, bam, bai -> [meta, bam, bai] }

    UMICOLLAPSE(bam_bai_ch, params.mode)
}
