#!/usr/bin/env nextflow
// Agent-facing wrapper adapter for local/modules/chai1.
// Follows docs/contracts/wrapper.md.
nextflow.enable.dsl = 2

include { CHAI1 } from '../main.nf'

params.inputs = null
params.sequences_tsv = null
params.cache_dir = null
params.target_sequence = null
params.ligand_smiles = null
params.use_msa_server = null
params.use_templates_server = null
params.num_trunk_recycles = null
params.num_diffn_samples = null
params.seed = null
params.args = null
params.outdir = null

workflow {
    def raw_inputs = params.inputs ? file(params.inputs, checkIfExists: true) : []
    def input_list = raw_inputs instanceof List ? raw_inputs : (raw_inputs ? [raw_inputs] : [])
    def seq_tsv = params.sequences_tsv ? file(params.sequences_tsv, checkIfExists: true) : []
    def cache = params.cache_dir ? file(params.cache_dir, checkIfExists: true) : []

    CHAI1(
        Channel.value([[id: 'chai1'], input_list, seq_tsv]),
        cache,
        params.target_sequence ?: '',
        params.ligand_smiles ?: '',
        params.use_msa_server,
        params.use_templates_server,
        params.num_trunk_recycles,
        params.num_diffn_samples,
        params.seed
    )
}
