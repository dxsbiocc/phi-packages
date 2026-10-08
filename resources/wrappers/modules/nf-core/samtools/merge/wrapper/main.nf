#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { SAMTOOLS_MERGE } from '../main.nf'

params.bam1   = null
params.bam2   = null
params.outdir = null

workflow {
    input_ch = Channel.value([
        [id: 'merged'],
        [file(params.bam1, checkIfExists: true), file(params.bam2, checkIfExists: true)],
        []
    ])
    fasta_ch = Channel.value([[:], [], [], []])

    SAMTOOLS_MERGE(input_ch, fasta_ch, '')
}
