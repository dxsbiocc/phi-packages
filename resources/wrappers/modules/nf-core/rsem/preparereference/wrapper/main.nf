#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf — see
// docs/design/phi-wrapper-agent-composition-design.md section 1. No
// `--star` in ext.args: builds RSEM's own (bowtie-based) reference index,
// not a combined STAR+RSEM one, matching the module's own default.
nextflow.enable.dsl = 2

include { RSEM_PREPAREREFERENCE } from '../main.nf'

params.fasta  = null
params.gtf    = null
params.outdir = null

workflow {
    RSEM_PREPAREREFERENCE(
        file(params.fasta, checkIfExists: true),
        file(params.gtf, checkIfExists: true)
    )
}
