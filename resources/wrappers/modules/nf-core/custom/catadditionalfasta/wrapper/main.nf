#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { CUSTOM_CATADDITIONALFASTA } from '../main.nf'

params.fasta      = null
params.gtf        = null
params.add_fasta  = null
params.biotype    = 'transgene'
params.outdir     = null

workflow {
    genome_ch = Channel
        .fromPath(params.fasta, checkIfExists: true)
        .map { fasta -> [[id: 'genome'], fasta] }
        .combine(Channel.fromPath(params.gtf, checkIfExists: true))
        .map { meta, fasta, gtf -> [meta, fasta, gtf] }

    add_fasta_ch = Channel
        .fromPath(params.add_fasta, checkIfExists: true)
        .map { fasta -> [[id: 'additional'], fasta] }

    CUSTOM_CATADDITIONALFASTA(genome_ch, add_fasta_ch, params.biotype)
}
