#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { BCFTOOLS_CSQ } from '../main.nf'


params.vcf        = null
params.fasta      = null
params.gff3       = null
params.outdir     = null

workflow {
    vcf_ch   = Channel.value([[id: 'result'], file(params.vcf, checkIfExists: true)])
    fasta_ch = Channel.value([[:], file(params.fasta, checkIfExists: true)])
    fai_ch   = Channel.value([[:], []])
    gff_ch   = Channel.value([[:], file(params.gff3, checkIfExists: true)])

    BCFTOOLS_CSQ(vcf_ch, fasta_ch, fai_ch, gff_ch)
}
