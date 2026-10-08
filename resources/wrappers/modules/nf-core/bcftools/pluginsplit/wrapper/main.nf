#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { BCFTOOLS_PLUGINSPLIT } from '../main.nf'


params.vcf        = null
params.tbi        = null
params.outdir     = null

workflow {
    ch_samples = Channel.of("normal\t-\tnormal", "tumour\t-\ttumour").collectFile(name: 'samples.txt', newLine: true)
    in_ch = Channel.of([[id: 'result'], file(params.vcf, checkIfExists: true), file(params.tbi, checkIfExists: true)])
                   .combine(ch_samples)
                   .combine(Channel.of([[], [], []]))

    BCFTOOLS_PLUGINSPLIT(in_ch)
}
