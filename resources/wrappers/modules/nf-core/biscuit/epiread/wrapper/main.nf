#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { BISCUIT_EPIREAD } from '../main.nf'
include { BISCUIT_INDEX } from '../../index/main.nf'

params.bam        = null
params.bai        = null
params.fasta      = null
params.outdir     = null

workflow {
    bam_ch = Channel.value([[id: 'sample'], file(params.bam, checkIfExists: true)])
    bai_ch = Channel.value([[id: 'sample'], file(params.bai, checkIfExists: true)])
    fasta_ch = Channel.value([[id: 'genome'], file(params.fasta, checkIfExists: true)])
    BISCUIT_INDEX(fasta_ch)
    snp_ch = Channel.value([[id: 'sample'], []])

    BISCUIT_EPIREAD(bam_ch, bai_ch, snp_ch, fasta_ch, BISCUIT_INDEX.out.index)
}
