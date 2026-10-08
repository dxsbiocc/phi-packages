#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { GATK4_SVCLUSTER } from '../main.nf'


params.vcf        = null
params.tbi        = null
params.ploidy     = null
params.fasta      = null
params.fai        = null
params.dict       = null
params.outdir     = null

workflow {
    vcf_ch = Channel.value([[id: 'clustered'], [file(params.vcf, checkIfExists: true)], [file(params.tbi, checkIfExists: true)]])
    ploidy_ch = Channel.value(file(params.ploidy, checkIfExists: true))
    fasta_ch = Channel.value(file(params.fasta, checkIfExists: true))
    fai_ch   = Channel.value(file(params.fai, checkIfExists: true))
    dict_ch  = Channel.value(file(params.dict, checkIfExists: true))

    GATK4_SVCLUSTER(vcf_ch, ploidy_ch, fasta_ch, fai_ch, dict_ch)
}
