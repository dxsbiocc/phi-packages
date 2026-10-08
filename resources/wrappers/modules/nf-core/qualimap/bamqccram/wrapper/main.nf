#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { QUALIMAP_BAMQCCRAM } from '../main.nf'

params.cram   = null
params.crai   = null
params.fasta  = null
params.fai    = null
params.outdir = null

workflow {
    cram_ch  = Channel.value([[id: 'test', single_end: false], file(params.cram, checkIfExists: true), file(params.crai, checkIfExists: true)])
    gff_ch   = Channel.value([])
    fasta_ch = Channel.value(file(params.fasta, checkIfExists: true))
    fai_ch   = Channel.value(file(params.fai, checkIfExists: true))

    QUALIMAP_BAMQCCRAM(cram_ch, gff_ch, fasta_ch, fai_ch)
}
