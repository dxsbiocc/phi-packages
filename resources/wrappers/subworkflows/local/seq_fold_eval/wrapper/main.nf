#!/usr/bin/env nextflow
// Agent-facing wrapper adapter for local/subworkflows/seq-fold-eval.
// Follows docs/contracts/wrapper.md.
nextflow.enable.dsl = 2

include { SEQ_FOLD_EVAL } from '../main.nf'

params.structures = null
params.predictor = null
params.mpnn_weights_dir = null
params.boltz_cache_dir = null
params.af3_model_dir = null
params.af3_db_dir = null
params.reference = null
params.mpnn_model_type = null
params.num_seq_per_pdb = null
params.mpnn_batch_size = null
params.mpnn_temperature = null
params.mpnn_seed = null
params.mpnn_output_mode = null
params.fixed_residues = null
params.chains_to_design = null
params.target_sequence = null
params.ligand_ccd = null
params.ligand_smiles = null
params.predict_affinity = null
params.use_msa_server = null
params.use_potentials = null
params.recycling_steps = null
params.diffusion_samples = null
params.boltz_output_format = null
params.af3_model_seeds = null
params.af3_run_data_pipeline = null
params.binder_chains = null
params.target_chains = null
params.contact_cutoff = null
params.hbond_cutoff = null
params.clash_cutoff = null
params.outdir = null

workflow {
    def meta = [id: 'seq_fold_eval']
    def struct_files = file(params.structures, checkIfExists: true)
    def struct_list = struct_files instanceof List ? struct_files : [struct_files]
    def mpnn_weights = params.mpnn_weights_dir ? file(params.mpnn_weights_dir, checkIfExists: true) : []
    def boltz_cache = params.boltz_cache_dir ? file(params.boltz_cache_dir, checkIfExists: true) : []
    def af3_model = params.af3_model_dir ? file(params.af3_model_dir, checkIfExists: true) : []
    def af3_db = params.af3_db_dir ? file(params.af3_db_dir, checkIfExists: true) : []
    def ref_file = params.reference ? file(params.reference, checkIfExists: true) : []

    if (params.predictor == 'alphafold3' && !params.af3_model_dir) {
        error "seq-fold-eval requires `af3_model_dir` when `predictor` is set to 'alphafold3'."
    }

    SEQ_FOLD_EVAL(
        Channel.value([meta, struct_list]),
        Channel.value(mpnn_weights),
        Channel.value(boltz_cache),
        Channel.value(af3_model),
        Channel.value(af3_db),
        Channel.value(ref_file),
        params.predictor ?: 'boltz2',
        params.mpnn_model_type,
        params.num_seq_per_pdb,
        params.mpnn_batch_size,
        params.mpnn_temperature,
        params.mpnn_seed,
        params.mpnn_output_mode,
        params.fixed_residues ?: '',
        params.chains_to_design ?: '',
        params.target_sequence ?: '',
        params.ligand_ccd ?: '',
        params.ligand_smiles ?: '',
        params.predict_affinity,
        params.use_msa_server,
        params.use_potentials,
        params.recycling_steps,
        params.diffusion_samples,
        params.boltz_output_format,
        params.af3_model_seeds ?: '1',
        params.af3_run_data_pipeline,
        params.binder_chains ?: '',
        params.target_chains ?: '',
        params.contact_cutoff,
        params.hbond_cutoff,
        params.clash_cutoff
    )
}
