#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored subworkflow at ../main.nf.
// Includes the subworkflow only — never another wrapper; the subworkflow
// itself references the already-migrated modules under resources/wrappers/modules.
// See docs/contracts/wrapper.md. Manually wired by scripts/wire-manual-subworkflows.py.
nextflow.enable.dsl = 2

include { ORFTABLE_FASTA_GTF_BUILDORFCATALOGUE } from '../main.nf'

params.orf_table = null
params.caller = "ribocode"
params.fasta = null
params.gtf = null
params.collapse = true
params.outdir = "results"

workflow {
    orf_ch = Channel.value([[id: 'orf'], file(params.orf_table, checkIfExists: true), params.caller])
    fasta_ch = Channel.value([[id: 'reference'], file(params.fasta, checkIfExists: true)])
    gtf_ch = Channel.value([[id: 'reference'], file(params.gtf, checkIfExists: true)])

    ORFTABLE_FASTA_GTF_BUILDORFCATALOGUE(orf_ch, fasta_ch, gtf_ch, params.collapse)
}
