#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { GATK4_LEARNREADORIENTATIONMODEL } from '../main.nf'


params.f1r2       = null
params.outdir     = null

workflow {
    f1r2_ch = Channel.value([[id: 'test'], [file(params.f1r2, checkIfExists: true)]])

    GATK4_LEARNREADORIENTATIONMODEL(f1r2_ch)
}
