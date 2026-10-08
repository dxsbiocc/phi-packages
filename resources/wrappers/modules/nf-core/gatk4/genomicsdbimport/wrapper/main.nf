#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { GATK4_GENOMICSDBIMPORT } from '../main.nf'


params.vcf        = null
params.tbi        = null
params.intervals  = null
params.outdir     = null

workflow {
    vcf_ch = Channel.value([[id: 'test'], file(params.vcf, checkIfExists: true), file(params.tbi, checkIfExists: true),
                            file(params.intervals, checkIfExists: true), [], []])

    GATK4_GENOMICSDBIMPORT(vcf_ch, false, false, false)
}
