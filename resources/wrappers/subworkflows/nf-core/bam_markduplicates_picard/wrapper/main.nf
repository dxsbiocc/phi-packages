#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored subworkflow at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1. Includes
// the subworkflow only — never another wrapper (section 6).
// The `.markdup` id suffix keeps Picard's output name different from the
// input's (the module errors out when they match).
nextflow.enable.dsl = 2

include { BAM_MARKDUPLICATES_PICARD } from '../main.nf'

params.bam    = null
params.outdir = null

workflow {
    bam_ch = Channel
        .fromPath(params.bam, checkIfExists: true)
        .map { bam -> [[id: "${bam.simpleName}.markdup"], bam] }

    fasta_fai_ch = Channel.value([[:], [], []])

    BAM_MARKDUPLICATES_PICARD(bam_ch, fasta_fai_ch)
}
