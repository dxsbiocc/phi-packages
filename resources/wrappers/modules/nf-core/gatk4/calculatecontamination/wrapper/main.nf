#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { GATK4_CALCULATECONTAMINATION } from '../main.nf'


params.pileup     = null
params.outdir     = null

workflow {
    pileup_ch = Channel.value([[id: 'test'], file(params.pileup, checkIfExists: true), []])

    GATK4_CALCULATECONTAMINATION(pileup_ch)
}
