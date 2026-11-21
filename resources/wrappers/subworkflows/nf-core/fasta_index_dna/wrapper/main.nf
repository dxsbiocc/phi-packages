#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored subworkflow at ../main.nf.
// Includes the subworkflow only — never another wrapper; the subworkflow
// itself references the already-migrated modules under resources/wrappers/modules.
// See docs/contracts/wrapper.md. Manually wired by scripts/wire-manual-subworkflows.py.
nextflow.enable.dsl = 2

include { FASTA_INDEX_DNA } from '../main.nf'

params.fasta = null
params.altliftover = null
params.aligner = "bwamem2"
params.outdir = "results"

workflow {
    fasta_ch = Channel.value([[id: 'reference'], file(params.fasta, checkIfExists: true)])
    alt_ch = params.altliftover ? Channel.value([[id: 'reference'], file(params.altliftover, checkIfExists: true)]) : Channel.empty()

    FASTA_INDEX_DNA(fasta_ch, alt_ch, params.aligner)
}
