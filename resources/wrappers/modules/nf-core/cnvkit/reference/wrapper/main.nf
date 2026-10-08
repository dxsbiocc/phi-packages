#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { CNVKIT_REFERENCE } from '../main.nf'


params.fasta      = null
params.targets    = null
params.antitargets = null
params.outdir     = null

workflow {
    CNVKIT_REFERENCE(file(params.fasta, checkIfExists: true), file(params.targets, checkIfExists: true), file(params.antitargets, checkIfExists: true))
}
