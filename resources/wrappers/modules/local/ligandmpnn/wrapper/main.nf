#!/usr/bin/env nextflow
// Agent-facing wrapper adapter for local/modules/ligandmpnn.
// Follows docs/contracts/wrapper.md.
nextflow.enable.dsl = 2

include { LIGANDMPNN } from '../main.nf'

params.structures = null
params.weights_dir = null
params.model_type = null
params.num_seq_per_pdb = null
params.batch_size = null
params.temperature = null
params.seed = null
params.output_mode = null
params.fixed_residues = null
params.chains_to_design = null
params.args = null
params.outdir = null

workflow {
    def struct_files = file(params.structures, checkIfExists: true)
    def struct_list = struct_files instanceof List ? struct_files : [struct_files]
    def weights = params.weights_dir ? file(params.weights_dir, checkIfExists: true) : []

    LIGANDMPNN(
        Channel.value([[id: 'ligandmpnn'], struct_list]),
        weights,
        params.model_type,
        params.num_seq_per_pdb,
        params.batch_size,
        params.temperature,
        params.seed,
        params.output_mode,
        params.fixed_residues ?: '',
        params.chains_to_design ?: ''
    )
}
