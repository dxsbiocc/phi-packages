#!/usr/bin/env nextflow
// Agent-facing wrapper adapter for local/modules/protein-metrics.
// Follows docs/contracts/wrapper.md.
nextflow.enable.dsl = 2

include { PROTEIN_METRICS } from '../main.nf'

params.structures = null
params.reference = null
params.binder_chains = null
params.target_chains = null
params.contact_cutoff = null
params.hbond_cutoff = null
params.clash_cutoff = null
params.args = null
params.outdir = null

workflow {
    def structure_files = file(params.structures, checkIfExists: true)
    def structure_list = structure_files instanceof List ? structure_files : [structure_files]

    PROTEIN_METRICS(
        Channel.value([[id: 'protein_metrics'], structure_list]),
        params.reference ? file(params.reference, checkIfExists: true) : [],
        params.binder_chains ?: '',
        params.target_chains ?: '',
        params.contact_cutoff,
        params.hbond_cutoff,
        params.clash_cutoff
    )
}
