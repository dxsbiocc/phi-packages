#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored subworkflow at ../main.nf.
// Includes the subworkflow only — never another wrapper; the subworkflow
// itself references the already-migrated modules under resources/wrappers/modules.
// See docs/contracts/wrapper.md. Manually wired by scripts/wire-manual-subworkflows.py.
nextflow.enable.dsl = 2

include { FASTA_CLEAN_FCS } from '../main.nf'

params.fasta = null
params.database = null
params.ramdisk_path = null
params.outdir = "results"

workflow {
    fasta_ch = Channel
        .fromPath(params.fasta, checkIfExists: true)
        .map { f -> [[id: f.simpleName], [f]] }
    db_ch = Channel.value(file(params.database, checkIfExists: true))

    FASTA_CLEAN_FCS(fasta_ch, db_ch, params.ramdisk_path ?: null)
}
