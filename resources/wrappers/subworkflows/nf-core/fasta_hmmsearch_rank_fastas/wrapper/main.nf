#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored subworkflow at ../main.nf.
// Includes the subworkflow only — never another wrapper; the subworkflow
// itself references the already-migrated modules under resources/wrappers/modules.
// See docs/contracts/wrapper.md. Manually wired by scripts/wire-manual-subworkflows.py.
nextflow.enable.dsl = 2

include { FASTA_HMMSEARCH_RANK_FASTAS } from '../main.nf'

params.hmms = null
params.fasta = null
params.save_domtblout = false
params.outdir = "results"

workflow {
    hmms_ch = Channel
        .fromPath(params.hmms, checkIfExists: true)
        .map { f -> [[id: f.simpleName], f] }
    fasta_ch = Channel.value(file(params.fasta, checkIfExists: true))

    FASTA_HMMSEARCH_RANK_FASTAS(hmms_ch, fasta_ch, params.save_domtblout)
}
