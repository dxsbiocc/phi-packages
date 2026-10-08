#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { VCF2MAF } from '../main.nf'


params.vcf        = null
params.fasta      = null
params.outdir     = null

workflow {
    in_ch = Channel.value([[id: 'sample'], file(params.vcf, checkIfExists: true)])

    VCF2MAF(in_ch, file(params.fasta, checkIfExists: true), [])
}
