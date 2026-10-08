#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { WHATSHAP_STATS } from '../main.nf'


params.vcf        = null
params.outdir     = null

workflow {
    vcf_ch = Channel.value([[id: 'result'], file(params.vcf, checkIfExists: true)])

    WHATSHAP_STATS(vcf_ch, true, true, true)
}
