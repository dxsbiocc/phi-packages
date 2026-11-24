#!/usr/bin/env nextflow
// Agent-facing wrapper adapter for local/modules/alphafold3.
// Follows docs/contracts/wrapper.md.
nextflow.enable.dsl = 2

include { ALPHAFOLD3 } from '../main.nf'

params.inputs = null
params.sequences_tsv = null
params.model_dir = null
params.db_dir = null
params.target_sequence = null
params.ligand_ccd = null
params.ligand_smiles = null
params.model_seeds = null
params.run_data_pipeline = null
params.num_diffusion_samples = null
params.args = null
params.outdir = null

workflow {
    def raw_inputs = params.inputs ? file(params.inputs, checkIfExists: true) : []
    def input_list = raw_inputs instanceof List ? raw_inputs : (raw_inputs ? [raw_inputs] : [])
    def seq_tsv = params.sequences_tsv ? file(params.sequences_tsv, checkIfExists: true) : []
    def model_dir = file(params.model_dir, checkIfExists: true)
    def db_dir = params.db_dir ? file(params.db_dir, checkIfExists: true) : []

    ALPHAFOLD3(
        Channel.value([[id: 'alphafold3'], input_list, seq_tsv]),
        model_dir,
        db_dir,
        params.target_sequence ?: '',
        params.ligand_ccd ?: '',
        params.ligand_smiles ?: '',
        params.model_seeds ?: '1',
        params.run_data_pipeline,
        params.num_diffusion_samples
    )
}
