#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored subworkflow at ../main.nf.
// Includes the subworkflow only — never another wrapper; the subworkflow
// itself references the already-migrated modules under resources/wrappers/modules.
// See docs/contracts/wrapper.md. Manually wired by scripts/wire-manual-subworkflows.py.
nextflow.enable.dsl = 2

include { FAA_SEQFU_SEQKIT } from '../main.nf'

params.fasta = null
params.skip_preprocessing = false
params.outdir = "results"

workflow {
    fasta_ch = Channel
        .fromPath(params.fasta, checkIfExists: true)
        .map { f -> [[id: f.simpleName], f] }

    FAA_SEQFU_SEQKIT(fasta_ch, params.skip_preprocessing)
}
