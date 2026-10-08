#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { SAMTOOLS_INDEX } from '../../../samtools/index/main.nf'
include { RSEQC_READDISTRIBUTION } from '../main.nf'

params.bam    = null
params.bed    = null
params.outdir = null

workflow {
    bam_ch = Channel
        .fromPath(params.bam, checkIfExists: true)
        .map { bam -> [[id: bam.simpleName], bam] }

    SAMTOOLS_INDEX(bam_ch)

    bam_bai_ch = bam_ch.join(SAMTOOLS_INDEX.out.index)

    bed_ch = Channel.fromPath(params.bed, checkIfExists: true)

    RSEQC_READDISTRIBUTION(bam_bai_ch, bed_ch)
}
