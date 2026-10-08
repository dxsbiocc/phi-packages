#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { SNPSIFT_ANNMEMCREATE } from '../main.nf'


params.vcf        = null
params.tbi        = null
params.fields     = 'DP,VDB'
params.outdir     = null

workflow {
    in_ch = Channel.value([[id: 'sample'], file(params.vcf, checkIfExists: true), file(params.tbi, checkIfExists: true), params.fields.tokenize(',')])

    SNPSIFT_ANNMEMCREATE(in_ch)
}
