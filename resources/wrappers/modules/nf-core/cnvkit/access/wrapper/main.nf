#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { CNVKIT_ACCESS } from '../main.nf'


params.fasta      = null
params.exclude    = null
params.outdir     = null

workflow {
    fasta_ch = Channel.value([[id: 'sample'], file(params.fasta, checkIfExists: true)])
    exclude_ch = Channel.value([[id: 'sample'], file(params.exclude, checkIfExists: true)])

    CNVKIT_ACCESS(fasta_ch, exclude_ch)
}
