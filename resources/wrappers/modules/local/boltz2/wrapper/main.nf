#!/usr/bin/env nextflow
// Agent-facing wrapper adapter for local/modules/boltz2.
// Follows docs/contracts/wrapper.md.
nextflow.enable.dsl = 2

include { BOLTZ2 } from '../main.nf'

params.inputs = null
params.sequences_tsv = null
params.cache_dir = null
params.target_sequence = null
params.ligand_ccd = null
params.ligand_smiles = null
params.predict_affinity = null
params.use_msa_server = null
params.use_potentials = null
params.recycling_steps = null
params.diffusion_samples = null
params.output_format = null
params.args = null
params.outdir = null

workflow {
    def raw_inputs = params.inputs ? file(params.inputs, checkIfExists: true) : []
    def input_list = raw_inputs instanceof List ? raw_inputs : (raw_inputs ? [raw_inputs] : [])
    def seq_tsv = params.sequences_tsv ? file(params.sequences_tsv, checkIfExists: true) : []
    def cache = params.cache_dir ? file(params.cache_dir, checkIfExists: true) : []

    BOLTZ2(
        Channel.value([[id: 'boltz2'], input_list, seq_tsv]),
        cache,
        params.target_sequence ?: '',
        params.ligand_ccd ?: '',
        params.ligand_smiles ?: '',
        params.predict_affinity,
        params.use_msa_server,
        params.use_potentials,
        params.recycling_steps,
        params.diffusion_samples,
        params.output_format
    )
}
