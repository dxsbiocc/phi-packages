#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { PICARD_COLLECTVARIANTCALLINGMETRICS } from '../main.nf'


params.vcf        = null
params.tbi        = null
params.fasta      = null
params.dict       = null
params.dbsnp      = null
params.dbsnp_tbi  = null
params.outdir     = null

workflow {
    input_ch = Channel.value([[id: 'test'], file(params.vcf, checkIfExists: true), file(params.tbi, checkIfExists: true), [],
                              file(params.fasta, checkIfExists: true), file(params.dict, checkIfExists: true),
                              file(params.dbsnp, checkIfExists: true), file(params.dbsnp_tbi, checkIfExists: true)])

    PICARD_COLLECTVARIANTCALLINGMETRICS(input_ch)
}
