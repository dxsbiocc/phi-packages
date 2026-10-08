#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { PRESEQ_LCEXTRAP } from '../main.nf'

params.bam        = null
params.single_end = true
params.outdir     = null

workflow {
    bam_ch = Channel
        .fromPath(params.bam, checkIfExists: true)
        .map { bam -> [[id: bam.simpleName, single_end: params.single_end], bam] }

    PRESEQ_LCEXTRAP(bam_ch)
}
