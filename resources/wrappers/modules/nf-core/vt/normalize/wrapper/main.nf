#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { VT_NORMALIZE } from '../main.nf'


params.vcf        = null
params.tbi        = null
params.fasta      = null
params.fai        = null
params.outdir     = null

workflow {
    vcf_ch   = Channel.value([[id: 'result', single_end: false], file(params.vcf, checkIfExists: true), file(params.tbi, checkIfExists: true), []])
    fasta_ch = Channel.value([[id: 'fasta'], file(params.fasta, checkIfExists: true)])
    fai_ch   = Channel.value([[id: 'fai'], file(params.fai, checkIfExists: true)])

    VT_NORMALIZE(vcf_ch, fasta_ch, fai_ch)
}
