#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { GATK4_UNMARKDUPLICATES } from '../main.nf'


params.bam        = null
params.outdir     = null

workflow {
    bam_ch = Channel.value([[id: 'unmarked', single_end: false], file(params.bam, checkIfExists: true)])

    GATK4_UNMARKDUPLICATES(bam_ch)
}
