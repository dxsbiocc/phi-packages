#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf — see
// docs/design/phi-wrapper-agent-composition-design.md section 1. No
// reference FASTA and no --write-index: plain coordinate-sort to BAM,
// matching the module's own defaults when those are left empty.
nextflow.enable.dsl = 2

include { SAMTOOLS_SORT } from '../main.nf'

params.bam    = null
params.outdir = null

workflow {
    Channel
        .fromPath(params.bam, checkIfExists: true)
        .map { bam -> [[id: bam.baseName], bam] }
        .set { bam_ch }

    ch_fasta = channel.value([[:], [], []])

    SAMTOOLS_SORT(bam_ch, ch_fasta, '')
}
