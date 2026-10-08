#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { GATK4_ASEREADCOUNTER } from '../main.nf'


params.bam        = null
params.bai        = null
params.vcf        = null
params.tbi        = null
params.fasta      = null
params.fai        = null
params.dict       = null
params.intervals  = null
params.outdir     = null

workflow {
    bam_ch = Channel.value([[id: 'test'], file(params.bam, checkIfExists: true), file(params.bai, checkIfExists: true),
                            file(params.vcf, checkIfExists: true), file(params.tbi, checkIfExists: true)])
    fasta_ch = Channel.value([[id: 'genome'], file(params.fasta, checkIfExists: true)])
    fai_ch   = Channel.value([[id: 'genome'], file(params.fai, checkIfExists: true)])
    dict_ch  = Channel.value([[id: 'genome'], file(params.dict, checkIfExists: true)])
    intervals_ch = Channel.value(file(params.intervals, checkIfExists: true))

    GATK4_ASEREADCOUNTER(bam_ch, fasta_ch, fai_ch, dict_ch, intervals_ch)
}
