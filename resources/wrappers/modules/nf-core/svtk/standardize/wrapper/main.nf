#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { SVTK_STANDARDIZE } from '../main.nf'


params.vcf        = null
params.fai        = null
params.outdir     = null

workflow {
    in_ch = Channel.value([[id: 'sample', caller: 'manta'], file(params.vcf, checkIfExists: true)])
    fai_ch = Channel.value([[id: 'fai'], file(params.fai, checkIfExists: true)])

    SVTK_STANDARDIZE(in_ch, fai_ch)
}
