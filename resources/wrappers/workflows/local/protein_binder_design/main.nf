//
// End-to-End Protein / Small-Molecule Binder Design Workflow
//
// Orchestrates:
//   1. RFD3 (parallelized across n_chunks with distinct chunk_id seeds)
//   2. SEQ_FOLD_EVAL (LIGANDMPNN -> BOLTZ2 / ALPHAFOLD3 / PROTENIX / CHAI1 -> PROTEIN_METRICS)
//   3. DESIGN_FILTER (campaign-wide lineage merge, boolean DSL filtering,
//      ranking, and per-parent diversity selection)
//

include { RFD3          } from '../../../modules/local/rfd3/main.nf'
include { SEQ_FOLD_EVAL } from '../../../subworkflows/local/seq_fold_eval/main.nf'
include { DESIGN_FILTER } from '../../../modules/local/design_filter/main.nf'

workflow PROTEIN_BINDER_DESIGN {
    take:
    ch_inputs_json           // channel: path(inputs_json) or []
    ch_input_structure       // channel: path(input_structure) or []
    ch_rfd3_weights          // channel: path(weights_dir) or []
    ch_mpnn_weights          // channel: path(weights_dir) or []
    ch_boltz_cache           // channel: path(cache_dir) or []
    ch_af3_model_dir         // channel: path(model_dir) or []
    ch_af3_db_dir            // channel: path(db_dir) or []
    ch_protenix_cache        // channel: path(cache_dir) or []
    ch_chai1_cache           // channel: path(cache_dir) or []
    ch_reference             // channel: path(reference) or []
    predictor                // val: string ('boltz2' | 'alphafold3' | 'protenix' | 'chai1')
    design_name              // val: string
    contig                   // val: string
    ligand                   // val: string
    n_chunks                 // val: int
    n_batches                // val: int
    diffusion_batch_size     // val: int
    num_timesteps            // val: int
    low_memory_mode          // val: boolean
    mpnn_model_type          // val: string
    num_seq_per_pdb          // val: int
    mpnn_batch_size          // val: int
    mpnn_temperature         // val: float
    mpnn_seed                // val: int
    mpnn_output_mode         // val: string
    fixed_residues           // val: string
    chains_to_design         // val: string
    target_sequence          // val: string
    ligand_ccd               // val: string
    ligand_smiles            // val: string
    predict_affinity         // val: boolean
    use_msa_server           // val: boolean
    use_potentials           // val: boolean
    recycling_steps          // val: int
    diffusion_samples        // val: int
    boltz_output_format      // val: string
    af3_model_seeds          // val: string
    af3_run_data_pipeline    // val: boolean
    protenix_model_name      // val: string
    protenix_diffusion_steps // val: int
    protenix_dtype           // val: string
    chai1_use_templates      // val: boolean
    binder_chains            // val: string
    target_chains            // val: string
    contact_cutoff           // val: float
    hbond_cutoff             // val: float
    clash_cutoff             // val: float
    filter_expr              // val: string
    rank_by                  // val: string
    rank_direction           // val: string
    group_by                 // val: string
    max_per_group            // val: int
    top_n                    // val: int

    main:
    ch_versions = Channel.empty()

    ch_rfd3_chunks = Channel
        .of(0 ..< (n_chunks as int))
        .combine(ch_inputs_json)
        .combine(ch_input_structure)
        .map { chunk_id, json_file, struct_file ->
            def meta = [id: "${design_name}_chunk${chunk_id}", chunk_id: chunk_id]
            [meta, json_file, struct_file, chunk_id]
        }

    RFD3(
        ch_rfd3_chunks.map { meta, json_file, struct_file, _cid -> [meta, json_file, struct_file] },
        ch_rfd3_weights,
        contig,
        ligand,
        n_batches,
        diffusion_batch_size,
        num_timesteps,
        low_memory_mode,
        ch_rfd3_chunks.map { _meta, _j, _s, cid -> cid }
    )
    ch_versions = ch_versions.mix(RFD3.out.versions)

    SEQ_FOLD_EVAL(
        RFD3.out.structures,
        ch_mpnn_weights,
        ch_boltz_cache,
        ch_af3_model_dir,
        ch_af3_db_dir,
        ch_protenix_cache,
        ch_chai1_cache,
        ch_reference,
        predictor,
        mpnn_model_type,
        num_seq_per_pdb,
        mpnn_batch_size,
        mpnn_temperature,
        mpnn_seed,
        mpnn_output_mode,
        fixed_residues,
        chains_to_design,
        target_sequence,
        ligand_ccd,
        ligand_smiles,
        predict_affinity,
        use_msa_server,
        use_potentials,
        recycling_steps,
        diffusion_samples,
        boltz_output_format,
        af3_model_seeds,
        af3_run_data_pipeline,
        protenix_model_name,
        protenix_diffusion_steps,
        protenix_dtype,
        chai1_use_templates,
        binder_chains,
        target_chains,
        contact_cutoff,
        hbond_cutoff,
        clash_cutoff
    )
    ch_versions = ch_versions.mix(SEQ_FOLD_EVAL.out.versions)

    ch_campaign_tables = SEQ_FOLD_EVAL.out.all_metric_tables
        .flatMap { _meta, table_list -> table_list }
        .collect()
        .map { all_tables -> [[id: design_name], all_tables] }

    DESIGN_FILTER(
        ch_campaign_tables,
        'design_id',
        'outer',
        filter_expr,
        rank_by,
        rank_direction,
        group_by,
        max_per_group,
        top_n
    )
    ch_versions = ch_versions.mix(DESIGN_FILTER.out.versions)

    emit:
    rfd3_structures      = RFD3.out.structures                   // channel: [ val(meta), [ path ] ]
    sequences_tsv        = SEQ_FOLD_EVAL.out.sequences_tsv       // channel: [ val(meta), path(tsv) ]
    predicted_structures = SEQ_FOLD_EVAL.out.predicted_structures// channel: [ val(meta), [ path ] ]
    merged_tsv           = DESIGN_FILTER.out.merged_tsv          // channel: [ val(meta), path(tsv) ]
    filtered_tsv         = DESIGN_FILTER.out.filtered_tsv        // channel: [ val(meta), path(tsv) ]
    summary_json         = DESIGN_FILTER.out.summary_json        // channel: [ val(meta), path(json) ]
    versions             = ch_versions                           // channel: [ path(versions.yml) ]
}
