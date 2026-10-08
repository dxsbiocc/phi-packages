#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { PREPROCESS_TRANSCRIPTS_FASTA_GENCODE } from '../main.nf'

params.fasta  = null
params.outdir = null

workflow {
    fasta_ch = Channel.fromPath(params.fasta, checkIfExists: true)

    PREPROCESS_TRANSCRIPTS_FASTA_GENCODE(fasta_ch)
}
