#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { GLIMPSE_CHUNK } from '../main.nf'


params.vcf        = null
params.csi        = null
params.region     = 'chr22'
params.outdir     = null

workflow {
    in_ch = Channel.value([[id: 'input'], file(params.vcf, checkIfExists: true), file(params.csi, checkIfExists: true), params.region])

    GLIMPSE_CHUNK(in_ch)
}
