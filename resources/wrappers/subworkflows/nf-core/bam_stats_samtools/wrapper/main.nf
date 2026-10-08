#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored subworkflow at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1. Includes
// the subworkflow only — never another wrapper (section 6).
nextflow.enable.dsl = 2

include { SAMTOOLS_INDEX } from '../../../../modules/nf-core/samtools/index/main'
include { BAM_STATS_SAMTOOLS } from '../main.nf'

params.bam        = null
params.single_end = false
params.outdir     = null

workflow {
    bam_ch = Channel
        .fromPath(params.bam, checkIfExists: true)
        .map { bam -> [[id: bam.simpleName, single_end: params.single_end], bam] }

    SAMTOOLS_INDEX(bam_ch)

    bam_bai_ch = bam_ch.join(SAMTOOLS_INDEX.out.index)

    fasta_fai_ch = Channel.value([[:], [], []])

    BAM_STATS_SAMTOOLS(bam_bai_ch, fasta_fai_ch)
}
