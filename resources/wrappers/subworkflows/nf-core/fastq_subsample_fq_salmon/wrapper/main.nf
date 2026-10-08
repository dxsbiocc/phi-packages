#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored subworkflow at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1. Includes
// the subworkflow only — never another wrapper (section 6).
nextflow.enable.dsl = 2

include { FASTQ_SUBSAMPLE_FQ_SALMON } from '../main.nf'

// Two explicit read-file params instead of one `reads` glob: Nextflow
// rejects glob patterns entirely for `https://` sources (same reasoning as
// the rsem/calculateexpression module wrapper).
params.reads_1             = null
params.reads_2             = null
params.genome_fasta        = null
params.transcript_fasta    = null
params.gtf                 = null
params.outdir              = null

workflow {
    reads_ch = Channel.of([
        [id: 'test', single_end: false],
        [file(params.reads_1, checkIfExists: true), file(params.reads_2, checkIfExists: true)]
    ])

    genome_fasta_ch     = Channel.of(file(params.genome_fasta, checkIfExists: true))
    transcript_fasta_ch = Channel.of(file(params.transcript_fasta, checkIfExists: true))
    gtf_ch              = Channel.of(file(params.gtf, checkIfExists: true))
    index_ch             = Channel.empty()

    FASTQ_SUBSAMPLE_FQ_SALMON(reads_ch, genome_fasta_ch, transcript_fasta_ch, gtf_ch, index_ch, true)
}
