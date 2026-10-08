#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { BCFTOOLS_MPILEUP } from '../main.nf'

params.bam         = null
params.fasta       = null
params.save_mpileup = false
params.outdir       = null

workflow {
    bam_ch = Channel
        .fromPath(params.bam, checkIfExists: true)
        .map { bam -> [[id: bam.simpleName], bam, [], []] }

    fasta_ch = Channel.value([[id: 'genome'], file(params.fasta, checkIfExists: true), []])

    BCFTOOLS_MPILEUP(bam_ch, fasta_ch, params.save_mpileup)
}
