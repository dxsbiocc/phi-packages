#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { GATK4_MODELSEGMENTS } from '../main.nf'
include { GATK4_PREPROCESSINTERVALS } from '../../preprocessintervals/main.nf'
include { GATK4_COLLECTREADCOUNTS } from '../../collectreadcounts/main.nf'
include { GATK4_CREATEREADCOUNTPANELOFNORMALS } from '../../createreadcountpanelofnormals/main.nf'
include { GATK4_DENOISEREADCOUNTS } from '../../denoisereadcounts/main.nf'

params.bam        = null
params.bai        = null
params.fasta      = null
params.fai        = null
params.dict       = null
params.outdir     = null

workflow {
    fasta_ch = Channel.value([[id: 'genome'], file(params.fasta, checkIfExists: true)])
    fai_ch   = Channel.value([[id: 'genome'], file(params.fai, checkIfExists: true)])
    dict_ch  = Channel.value([[id: 'genome'], file(params.dict, checkIfExists: true)])
    empty_ch = Channel.value([[:], []])
    GATK4_PREPROCESSINTERVALS(fasta_ch, fai_ch, dict_ch, empty_ch, empty_ch)
    bam_ch = Channel.value([[id: 'test', single_end: false], file(params.bam, checkIfExists: true), file(params.bai, checkIfExists: true)])
    GATK4_COLLECTREADCOUNTS(bam_ch.combine(GATK4_PREPROCESSINTERVALS.out.interval_list.map { meta, list -> list }), fasta_ch, fai_ch, dict_ch)
    GATK4_CREATEREADCOUNTPANELOFNORMALS(GATK4_COLLECTREADCOUNTS.out.tsv.groupTuple())
    GATK4_DENOISEREADCOUNTS(GATK4_COLLECTREADCOUNTS.out.tsv.first(), GATK4_CREATEREADCOUNTPANELOFNORMALS.out.pon)
    GATK4_MODELSEGMENTS(GATK4_DENOISEREADCOUNTS.out.denoised.first())
}
