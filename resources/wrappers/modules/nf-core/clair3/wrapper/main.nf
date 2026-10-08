#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { CLAIR3 } from '../main.nf'
include { UNTAR } from '../../untar/main.nf'

params.bam        = null
params.bai        = null
params.model      = null
params.fasta      = null
params.fai        = null
params.outdir     = null

workflow {
    UNTAR(Channel.value([[id: 'model'], file(params.model, checkIfExists: true)]))
    in_ch = Channel.value([[id: 'test'], file(params.bam, checkIfExists: true), file(params.bai, checkIfExists: true), []])
        .join(UNTAR.out.untar)
        .combine(Channel.of('hifi'))
    fasta_ch = Channel.value([[id: 'test'], file(params.fasta, checkIfExists: true)])
    fai_ch = Channel.value([[id: 'test'], file(params.fai, checkIfExists: true)])

    CLAIR3(in_ch, fasta_ch, fai_ch)
}
