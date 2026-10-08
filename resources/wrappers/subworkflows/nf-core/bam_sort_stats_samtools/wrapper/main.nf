#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored subworkflow at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1. Includes
// the subworkflow only — never another wrapper (section 6).
nextflow.enable.dsl = 2

include { BAM_SORT_STATS_SAMTOOLS } from '../main.nf'

params.bam    = null
params.outdir = null

workflow {
    // The `.sorted` id suffix keeps SAMTOOLS_SORT's output name different from the input's.
    bam_ch = Channel
        .fromPath(params.bam, checkIfExists: true)
        .map { bam -> [[id: "${bam.baseName}.sorted", single_end: false], bam] }

    ch_fasta_fai = channel.value([[:], [], []])

    BAM_SORT_STATS_SAMTOOLS(bam_ch, ch_fasta_fai)
}
