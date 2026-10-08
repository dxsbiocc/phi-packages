#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { BCFTOOLS_CONSENSUS } from '../main.nf'


params.vcf        = null
params.tbi        = null
params.fasta      = null
params.mask       = null
params.outdir     = null

workflow {
    in_ch = Channel.value([[id: 'result'], file(params.vcf, checkIfExists: true), file(params.tbi, checkIfExists: true),
                           file(params.fasta, checkIfExists: true), file(params.mask, checkIfExists: true)])

    BCFTOOLS_CONSENSUS(in_ch)
}
