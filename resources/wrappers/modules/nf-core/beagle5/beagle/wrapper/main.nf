#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { BEAGLE5_BEAGLE } from '../main.nf'


params.vcf        = null
params.csi        = null
params.outdir     = null

workflow {
    in_ch = Channel.value([[id: 'sample'], file(params.vcf, checkIfExists: true), file(params.csi, checkIfExists: true), [], [], [], [], [], []])

    BEAGLE5_BEAGLE(in_ch)
}
