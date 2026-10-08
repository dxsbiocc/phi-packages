#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { GATK_INDELREALIGNER } from '../main.nf'
include { GATK_REALIGNERTARGETCREATOR } from '../../realignertargetcreator/main.nf'

params.bam        = null
params.bai        = null
params.fasta      = null
params.fai        = null
params.dict       = null
params.outdir     = null

workflow {
    bam_ch = Channel.value([[id: 'test'], file(params.bam, checkIfExists: true), file(params.bai, checkIfExists: true)])
    fasta_ch = Channel.value([[id: 'test'], file(params.fasta, checkIfExists: true)])
    fai_ch   = Channel.value([[id: 'test'], file(params.fai, checkIfExists: true)])
    dict_ch  = Channel.value([[id: 'test'], file(params.dict, checkIfExists: true)])

    GATK_REALIGNERTARGETCREATOR(bam_ch, fasta_ch, fai_ch, dict_ch, [[], []])
    GATK_INDELREALIGNER(bam_ch.join(GATK_REALIGNERTARGETCREATOR.out.intervals), fasta_ch, fai_ch, dict_ch, [[], []])
}
