#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// Composes the vendored samtools/fixmate and samtools/sort modules
// (../../fixmate/main.nf, ../../sort/main.nf) so the agent can mark
// duplicates on a raw BAM in one step — samtools markdup requires MC
// tags from fixmate and coordinate-sorted input. See
// docs/design/phi-wrapper-agent-composition-design.md sections 1 and 3.
nextflow.enable.dsl = 2

include { SAMTOOLS_FIXMATE } from '../../fixmate/main.nf'
include { SAMTOOLS_SORT }    from '../../sort/main.nf'
include { SAMTOOLS_MARKDUP } from '../main.nf'

params.bam    = null
params.outdir = null

workflow {
    bam_ch = Channel.value([[id: 'test'], file(params.bam, checkIfExists: true)])
    fasta_ch = Channel.value([[:], [], []])

    SAMTOOLS_FIXMATE(bam_ch, fasta_ch)
    SAMTOOLS_SORT(SAMTOOLS_FIXMATE.out.bam, fasta_ch, '')
    SAMTOOLS_MARKDUP(SAMTOOLS_SORT.out.bam, fasta_ch)
}
