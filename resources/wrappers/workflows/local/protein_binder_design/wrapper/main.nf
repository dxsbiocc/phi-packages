#!/usr/bin/env nextflow
// Agent-facing wrapper adapter for local/workflows/protein-binder-design.
// Follows docs/contracts/wrapper.md.
nextflow.enable.dsl = 2

include { PROTEIN_BINDER_DESIGN } from '../main.nf'

params.design_name = null
params.predictor = null
params.inputs_json = null
params.input_structure = null
params.contig = null
params.ligand = null
params.rfd3_weights_dir = null
params.mpnn_weights_dir = null
params.boltz_cache_dir = null
params.af3_model_dir = null
params.af3_db_dir = null
params.protenix_cache_dir = null
params.chai1_cache_dir = null
params.reference = null
params.n_chunks = null
params.n_batches = null
params.diffusion_batch_size = null
params.num_timesteps = null
params.low_memory_mode = null
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
params.protenix_model_name = null
params.protenix_diffusion_steps = null
params.protenix_dtype = null
params.chai1_use_templates = null
params.binder_chains = null
params.target_chains = null
params.contact_cutoff = null
params.hbond_cutoff = null
params.clash_cutoff = null
params.run_foldseek_cluster = null
params.foldseek_coverage = null
params.foldseek_tmscore_threshold = null
params.filter_expr = null
params.rank_by = null
params.rank_direction = null
params.group_by = null
params.max_per_group = null
params.top_n = null
params.outdir = null

workflow {
    if (params.predictor == 'alphafold3' && !params.af3_model_dir) {
        error "protein-binder-design requires `af3_model_dir` when `predictor` is set to 'alphafold3'."
    }

    def json_input = params.inputs_json ? file(params.inputs_json, checkIfExists: true) : []
    def struct_input = params.input_structure ? file(params.input_structure, checkIfExists: true) : []
    def rfd3_weights = params.rfd3_weights_dir ? file(params.rfd3_weights_dir, checkIfExists: true) : []
    def mpnn_weights = params.mpnn_weights_dir ? file(params.mpnn_weights_dir, checkIfExists: true) : []
    def boltz_cache = params.boltz_cache_dir ? file(params.boltz_cache_dir, checkIfExists: true) : []
    def af3_model = params.af3_model_dir ? file(params.af3_model_dir, checkIfExists: true) : []
    def af3_db = params.af3_db_dir ? file(params.af3_db_dir, checkIfExists: true) : []
    def protenix_cache = params.protenix_cache_dir ? file(params.protenix_cache_dir, checkIfExists: true) : []
    def chai1_cache = params.chai1_cache_dir ? file(params.chai1_cache_dir, checkIfExists: true) : []
    def ref_file = params.reference ? file(params.reference, checkIfExists: true) : []

    PROTEIN_BINDER_DESIGN(
        Channel.value(json_input),
        Channel.value(struct_input),
        Channel.value(rfd3_weights),
        Channel.value(mpnn_weights),
        Channel.value(boltz_cache),
        Channel.value(af3_model),
        Channel.value(af3_db),
        Channel.value(protenix_cache),
        Channel.value(chai1_cache),
        Channel.value(ref_file),
        params.predictor ?: 'boltz2',
        params.design_name ?: 'binder_campaign',
        params.contig ?: '',
        params.ligand ?: '',
        params.n_chunks,
        params.n_batches,
        params.diffusion_batch_size,
        params.num_timesteps,
        params.low_memory_mode,
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
        params.protenix_model_name ?: 'protenix_base_default_v1.0.0',
        params.protenix_diffusion_steps,
        params.protenix_dtype ?: 'bf16',
        params.chai1_use_templates,
        params.binder_chains ?: '',
        params.target_chains ?: '',
        params.contact_cutoff,
        params.hbond_cutoff,
        params.clash_cutoff,
        params.run_foldseek_cluster,
        params.foldseek_coverage,
        params.foldseek_tmscore_threshold,
        params.filter_expr ?: '',
        params.rank_by ?: '',
        params.rank_direction,
        params.group_by ?: '',
        params.max_per_group,
        params.top_n
    )
}
