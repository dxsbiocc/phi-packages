#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored subworkflow at ../main.nf.
// Includes the subworkflow only — never another wrapper; the subworkflow
// itself references the already-migrated modules under resources/wrappers/modules.
// See docs/contracts/wrapper.md. Manually wired by scripts/wire-manual-subworkflows.py.
nextflow.enable.dsl = 2

include { MMSEQS_FASTA_CLUSTER } from '../main.nf'

params.fasta = null
params.clustering_tool = "linclust"
params.outdir = "results"

workflow {
    seq_ch = Channel
        .fromPath(params.fasta, checkIfExists: true)
        .map { f -> [[id: f.simpleName], f] }

    MMSEQS_FASTA_CLUSTER(seq_ch, params.clustering_tool)
}
