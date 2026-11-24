#!/usr/bin/env nextflow
// Agent-facing wrapper adapter for local/modules/protenix.
// Follows docs/contracts/wrapper.md.
nextflow.enable.dsl = 2

include { PROTENIX } from '../main.nf'

params.inputs = null
params.sequences_tsv = null
params.cache_dir = null
params.model_name = null
params.target_sequence = null
params.ligand_ccd = null
params.ligand_smiles = null
params.seeds = null
params.recycling_steps = null
params.diffusion_steps = null
params.diffusion_samples = null
params.use_msa = null
params.dtype = null
params.args = null
params.outdir = null

workflow {
    def raw_inputs = params.inputs ? file(params.inputs, checkIfExists: true) : []
    def input_list = raw_inputs instanceof List ? raw_inputs : (raw_inputs ? [raw_inputs] : [])
    def seq_tsv = params.sequences_tsv ? file(params.sequences_tsv, checkIfExists: true) : []
    def cache = params.cache_dir ? file(params.cache_dir, checkIfExists: true) : []

    PROTENIX(
        Channel.value([[id: 'protenix'], input_list, seq_tsv]),
        cache,
        params.model_name ?: 'protenix_base_default_v1.0.0',
        params.target_sequence ?: '',
        params.ligand_ccd ?: '',
        params.ligand_smiles ?: '',
        params.seeds ?: '101',
        params.recycling_steps,
        params.diffusion_steps,
        params.diffusion_samples,
        params.use_msa,
        params.dtype ?: 'bf16'
    )
}
