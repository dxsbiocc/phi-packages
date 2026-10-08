#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { VARSCAN_PROCESSSOMATIC } from '../main.nf'
include { VARSCAN_SOMATIC } from '../../somatic/main.nf'

params.normal     = null
params.tumour     = null
params.outdir     = null

workflow {
    in_ch = Channel.value([[id: 'sample', single_end: false], file(params.normal, checkIfExists: true), file(params.tumour, checkIfExists: true)])

    VARSCAN_SOMATIC(in_ch)
    VARSCAN_PROCESSSOMATIC(VARSCAN_SOMATIC.out.vcf_snvs.collect { meta, vcf -> [meta, vcf] })
}
