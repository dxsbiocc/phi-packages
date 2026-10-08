#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { DUPRADAR } from '../main.nf'

params.bam          = null
params.gtf          = null
params.single_end   = true
params.strandedness = 'forward'
params.outdir       = null

workflow {
    bam_ch = Channel
        .fromPath(params.bam, checkIfExists: true)
        .map { bam -> [[id: bam.simpleName, single_end: params.single_end, strandedness: params.strandedness], bam] }

    gtf_ch = Channel
        .fromPath(params.gtf, checkIfExists: true)
        .map { gtf -> [[id: gtf.simpleName], gtf] }

    DUPRADAR(bam_ch, gtf_ch)
}
