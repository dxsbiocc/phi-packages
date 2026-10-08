#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { QUALIMAP_RNASEQ } from '../main.nf'

params.bam          = null
params.gtf          = null
params.single_end   = false
params.strandedness = 'unstranded'
params.outdir       = null

workflow {
    bam_ch = Channel
        .fromPath(params.bam, checkIfExists: true)
        .map { bam -> [[id: bam.simpleName, single_end: params.single_end, strandedness: params.strandedness], bam] }

    gtf_ch = Channel
        .fromPath(params.gtf, checkIfExists: true)
        .map { gtf -> [[id: 'gtf'], gtf] }

    QUALIMAP_RNASEQ(bam_ch, gtf_ch)
}
