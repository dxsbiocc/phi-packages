#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
// Manually wired: the module's staging pattern is too exotic for the generator.
nextflow.enable.dsl = 2

include { METABULI_BUILD } from '../main.nf'

params.fasta = null
params.taxonomy_names = null
params.taxonomy_nodes = null
params.taxonomy_merged = null
params.accession2taxid = null
params.cds_info = null
params.args = null
params.outdir = null

workflow {
    def opt = { p -> p ? file(p, checkIfExists: true) : [] }

    METABULI_BUILD(
        Channel.value([[id: 'metabuli'], [file(params.fasta, checkIfExists: true)]]),
        file(params.taxonomy_names, checkIfExists: true),
        file(params.taxonomy_nodes, checkIfExists: true),
        opt(params.taxonomy_merged),
        [file(params.accession2taxid, checkIfExists: true)],
        opt(params.cds_info)
    )
}
