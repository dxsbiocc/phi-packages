#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored subworkflow at ../main.nf.
// Includes the subworkflow only — never another wrapper; the subworkflow
// itself references the already-migrated modules under resources/wrappers/modules.
// See docs/contracts/wrapper.md. Manually wired by scripts/wire-manual-subworkflows.py.
nextflow.enable.dsl = 2

include { HOMER_GROSEQ } from '../main.nf'

params.bam = null
params.fasta = null
params.uniqmap = null
params.outdir = "results"

workflow {
    bam_ch = Channel
        .fromPath(params.bam, checkIfExists: true)
        .map { f -> [[id: f.simpleName], [f]] }
    fasta_ch = Channel.value(file(params.fasta, checkIfExists: true))

    HOMER_GROSEQ(bam_ch, fasta_ch, params.uniqmap ?: null)
}
