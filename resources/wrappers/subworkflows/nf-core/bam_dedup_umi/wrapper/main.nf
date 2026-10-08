#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored subworkflow at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1. Includes
// the subworkflow only — never another wrapper (section 6).
//
// Genome and transcriptome BAM both default to the same input file — as
// the subworkflow's own test does — since this wrapper's job is to exercise
// the composition, not model a real dual-alignment output.
nextflow.enable.dsl = 2

include { BAM_DEDUP_UMI } from '../main.nf'

params.genome_bam                  = null
params.genome_bai                  = null
params.fasta                       = null
params.fai                         = null
params.transcriptome_bam           = null
params.transcript_fasta            = null
params.umi_dedup_tool              = 'umitools'
params.umitools_dedup_stats        = false
params.umitools_dedup_primary_only = false
params.outdir                      = null

workflow {
    genome_bam_ch = Channel.value([
        [id: 'test'],
        file(params.genome_bam, checkIfExists: true),
        file(params.genome_bai, checkIfExists: true)
    ])

    fasta_fai_ch = Channel.value([
        [id: 'genome'],
        file(params.fasta, checkIfExists: true),
        file(params.fai, checkIfExists: true)
    ])

    transcriptome_bam_ch = Channel.value([
        [id: 'test'],
        file(params.transcriptome_bam, checkIfExists: true)
    ])

    transcript_fasta_fai_ch = Channel.value([
        [id: 'genome'],
        file(params.transcript_fasta, checkIfExists: true),
        []
    ])

    BAM_DEDUP_UMI(
        genome_bam_ch,
        fasta_fai_ch,
        params.umi_dedup_tool,
        params.umitools_dedup_stats,
        transcriptome_bam_ch,
        transcript_fasta_fai_ch,
        params.umitools_dedup_primary_only
    )
}
