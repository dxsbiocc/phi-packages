#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored subworkflow at ../main.nf.
// Includes the subworkflow only — never another wrapper; the subworkflow
// itself references the already-migrated modules under resources/wrappers/modules.
// See docs/contracts/wrapper.md. Manually wired by scripts/wire-manual-subworkflows.py.
nextflow.enable.dsl = 2

include { FASTA_CONSENSUS_AUTOCYCLER } from '../main.nf'

params.contigs = null
params.outdir = "results"

workflow {
    contigs_ch = Channel
        .fromPath(params.contigs, checkIfExists: true)
        .collect()
        .map { files -> [[id: 'consensus'], files] }

    FASTA_CONSENSUS_AUTOCYCLER(contigs_ch)
}
