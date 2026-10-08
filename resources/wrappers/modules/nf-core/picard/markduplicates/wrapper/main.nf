#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
// The `.markdup` id suffix keeps the output name different from the input's
// (the module errors out when they match).
nextflow.enable.dsl = 2

include { PICARD_MARKDUPLICATES } from '../main.nf'

params.bam    = null
params.outdir = null

workflow {
    bam_ch = Channel
        .fromPath(params.bam, checkIfExists: true)
        .map { bam -> [[id: "${bam.baseName}.markdup"], bam] }

    ch_fasta = channel.value([[:], [], []])

    PICARD_MARKDUPLICATES(bam_ch, ch_fasta)
}
