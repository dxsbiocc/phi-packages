#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { GATK4_ADDORREPLACEREADGROUPS } from '../main.nf'


params.bam        = null
params.outdir     = null

workflow {
    bam_ch = Channel.value([[id: 'replaced'], file(params.bam, checkIfExists: true)])
    ref_ch = Channel.value([[:], []])

    GATK4_ADDORREPLACEREADGROUPS(bam_ch, ref_ch, ref_ch)
}
