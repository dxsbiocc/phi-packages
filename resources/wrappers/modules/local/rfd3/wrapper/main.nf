#!/usr/bin/env nextflow
// Agent-facing wrapper adapter for local/modules/rfd3.
// Follows docs/contracts/wrapper.md.
nextflow.enable.dsl = 2

include { RFD3 } from '../main.nf'

params.inputs_json = null
params.input_structure = null
params.weights_dir = null
params.contig = null
params.ligand = null
params.n_batches = null
params.diffusion_batch_size = null
params.num_timesteps = null
params.low_memory_mode = null
params.chunk_id = null
params.args = null
params.outdir = null

workflow {
    def json_input = params.inputs_json ? file(params.inputs_json, checkIfExists: true) : []
    def struct_input = params.input_structure ? file(params.input_structure, checkIfExists: true) : []
    def weights = params.weights_dir ? file(params.weights_dir, checkIfExists: true) : []

    RFD3(
        Channel.value([[id: 'rfd3'], json_input, struct_input]),
        weights,
        params.contig ?: '',
        params.ligand ?: '',
        params.n_batches,
        params.diffusion_batch_size,
        params.num_timesteps,
        params.low_memory_mode,
        params.chunk_id
    )
}
