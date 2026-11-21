#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored subworkflow at ../main.nf.
// Includes the subworkflow only — never another wrapper; the subworkflow
// itself references the already-migrated modules under resources/wrappers/modules.
// See docs/contracts/wrapper.md. Manually wired by scripts/wire-manual-subworkflows.py.
nextflow.enable.dsl = 2

include { FASTA_CLEAN_FAIDX } from '../main.nf'

params.fasta = null
params.replace_dots = false
params.get_chromsizes = true
params.get_dict = true
params.outdir = "results"

workflow {
    reference_ch = Channel
        .fromPath(params.fasta, checkIfExists: true)
        .map { f -> [[id: f.simpleName], f] }

    FASTA_CLEAN_FAIDX(reference_ch, params.replace_dots, params.get_chromsizes, params.get_dict)
}
