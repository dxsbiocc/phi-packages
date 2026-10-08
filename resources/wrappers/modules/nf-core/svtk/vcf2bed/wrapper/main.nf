#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { SVTK_VCF2BED } from '../main.nf'


params.vcf        = null
params.tbi        = null
params.outdir     = null

workflow {
    in_ch = Channel.value([[id: 'sample', single_end: false], file(params.vcf, checkIfExists: true), file(params.tbi, checkIfExists: true)])

    SVTK_VCF2BED(in_ch)
}
