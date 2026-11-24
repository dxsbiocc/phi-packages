#!/usr/bin/env nextflow
// Agent-facing wrapper adapter for local/modules/bindcraft.
// Follows docs/contracts/wrapper.md.
nextflow.enable.dsl = 2

include { BINDCRAFT } from '../main.nf'

params.starting_pdb = null
params.settings_json = null
params.filters_json = null
params.advanced_json = null
params.binder_name = null
params.chains = null
params.target_hotspot_residues = null
params.min_length = null
params.max_length = null
params.number_of_final_designs = null
params.args = null
params.outdir = null

workflow {
    def pdb = params.starting_pdb ? file(params.starting_pdb, checkIfExists: true) : []
    def settings = params.settings_json ? file(params.settings_json, checkIfExists: true) : []
    def filters = params.filters_json ? file(params.filters_json, checkIfExists: true) : []
    def advanced = params.advanced_json ? file(params.advanced_json, checkIfExists: true) : []

    BINDCRAFT(
        Channel.value([[id: params.binder_name ?: 'bindcraft'], pdb, settings, filters, advanced]),
        params.binder_name ?: 'binder',
        params.chains ?: 'A',
        params.target_hotspot_residues ?: '',
        params.min_length,
        params.max_length,
        params.number_of_final_designs
    )
}
