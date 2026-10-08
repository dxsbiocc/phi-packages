#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { GATK4_SITEDEPTHTOBAF } from '../main.nf'
include { GATK4_COLLECTSVEVIDENCE } from '../../collectsvevidence/main.nf'

params.bam        = null
params.bai        = null
params.bam2       = null
params.bai2       = null
params.vcf        = null
params.tbi        = null
params.fasta      = null
params.fai        = null
params.dict       = null
params.outdir     = null

workflow {
    fasta_ch = Channel.value(file(params.fasta, checkIfExists: true))
    fai_ch   = Channel.value(file(params.fai, checkIfExists: true))
    dict_ch  = Channel.value(file(params.dict, checkIfExists: true))
    vcf = file(params.vcf, checkIfExists: true)
    tbi = file(params.tbi, checkIfExists: true)
    GATK4_COLLECTSVEVIDENCE(Channel.of(
        [[id: 'tumor', single_end: false], file(params.bam, checkIfExists: true), file(params.bai, checkIfExists: true), vcf, tbi],
        [[id: 'normal', single_end: false], file(params.bam2, checkIfExists: true), file(params.bai2, checkIfExists: true), vcf, tbi]),
        fasta_ch, fai_ch, dict_ch)
    sd_ch = GATK4_COLLECTSVEVIDENCE.out.site_depths.combine(GATK4_COLLECTSVEVIDENCE.out.site_depths_index, by: 0)
                .map { meta, file, index -> [[id: 'test'], file, index] }.groupTuple()
    GATK4_SITEDEPTHTOBAF(sd_ch, [vcf, tbi], fasta_ch, fai_ch, dict_ch)
}
