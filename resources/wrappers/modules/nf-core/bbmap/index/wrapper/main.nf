#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { BBMAP_INDEX } from '../main.nf'


params.fasta      = null
params.outdir     = null

workflow {
    BBMAP_INDEX(file(params.fasta, checkIfExists: true))
}
