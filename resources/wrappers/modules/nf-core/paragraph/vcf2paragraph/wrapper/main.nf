#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { PARAGRAPH_VCF2PARAGRAPH } from '../main.nf'


params.vcf        = null
params.fasta      = null
params.outdir     = null

workflow {
    in_ch = Channel.value([[id: 'sample', single_end: false], file(params.vcf, checkIfExists: true)])
    fasta_ch = Channel.value([[], file(params.fasta, checkIfExists: true)])

    PARAGRAPH_VCF2PARAGRAPH(in_ch, fasta_ch)
}
