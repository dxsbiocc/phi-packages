#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { PLINK2_VCF2BGEN } from '../main.nf'


params.vcf        = null
params.outdir     = null

workflow {
    in_ch = Channel.value([[id: 'test'], file(params.vcf, checkIfExists: true), 'DS', true, 'double-id'])
    PLINK2_VCF2BGEN(in_ch)

}
