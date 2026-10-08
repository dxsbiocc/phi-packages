#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { SAMTOOLS_INDEX } from '../../samtools/index/main.nf'
include { RUSTQC } from '../main.nf'

params.bam          = null
params.gtf          = null
params.single_end   = false
params.outdir       = null

workflow {
    bam_ch = Channel
        .fromPath(params.bam, checkIfExists: true)
        .map { bam -> [[id: bam.simpleName, single_end: params.single_end], bam] }

    SAMTOOLS_INDEX(bam_ch)

    bam_bai_ch = bam_ch.join(SAMTOOLS_INDEX.out.index)

    gtf_ch = Channel
        .fromPath(params.gtf, checkIfExists: true)
        .map { gtf -> [[id: 'gtf'], gtf] }

    RUSTQC(bam_bai_ch, gtf_ch)
}
