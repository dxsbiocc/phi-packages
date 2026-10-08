#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf — see
// docs/design/phi-wrapper-agent-composition-design.md section 1. No index
// and no reference FASTA: `samtools stats` reads the BAM sequentially for
// whole-file stats without needing either, matching the module's defaults
// when those are left empty.
nextflow.enable.dsl = 2

include { SAMTOOLS_STATS } from '../main.nf'

params.bam    = null
params.outdir = null

workflow {
    Channel
        .fromPath(params.bam, checkIfExists: true)
        .map { bam -> [[id: bam.baseName], bam, []] }
        .set { bam_ch }

    ch_fasta = channel.value([[:], [], []])

    SAMTOOLS_STATS(bam_ch, ch_fasta)
}
