//
// Sequence Design -> Structure Prediction -> Structural Evaluation Subworkflow
//
// Chains LIGANDMPNN (inverse folding) -> selectable structure prediction oracle
// (BOLTZ2, ALPHAFOLD3, PROTENIX, or CHAI1 via `predictor`) -> PROTEIN_METRICS
// while preserving design_id / parent_design_id lineage across all stages.
//

include { LIGANDMPNN      } from '../../../modules/local/ligandmpnn/main.nf'
include { BOLTZ2          } from '../../../modules/local/boltz2/main.nf'
include { ALPHAFOLD3      } from '../../../modules/local/alphafold3/main.nf'
include { PROTENIX        } from '../../../modules/local/protenix/main.nf'
include { CHAI1           } from '../../../modules/local/chai1/main.nf'
include { PROTEIN_METRICS } from '../../../modules/local/protein_metrics/main.nf'

workflow SEQ_FOLD_EVAL {
    take:
    ch_structures          // channel: [ val(meta), [ path(structure), ... ] ]
    ch_mpnn_weights        // channel: path(weights_dir) or []
    ch_boltz_cache         // channel: path(cache_dir) or []
    ch_af3_model_dir       // channel: path(model_dir) or []
    ch_af3_db_dir          // channel: path(db_dir) or []
    ch_protenix_cache      // channel: path(cache_dir) or []
    ch_chai1_cache         // channel: path(cache_dir) or []
    ch_reference           // channel: path(reference) or []
    predictor              // val: string ('boltz2' | 'alphafold3' | 'protenix' | 'chai1')
    mpnn_model_type        // val: string
    num_seq_per_pdb        // val: int
    mpnn_batch_size        // val: int
    mpnn_temperature       // val: float
    mpnn_seed              // val: int
    mpnn_output_mode       // val: string
    fixed_residues         // val: string
    chains_to_design       // val: string
    target_sequence        // val: string
    ligand_ccd             // val: string
    ligand_smiles          // val: string
    predict_affinity       // val: boolean
    use_msa_server         // val: boolean
    use_potentials         // val: boolean
    recycling_steps        // val: int
    diffusion_samples      // val: int
    boltz_output_format    // val: string
    af3_model_seeds        // val: string
    af3_run_data_pipeline  // val: boolean
    protenix_model_name    // val: string
    protenix_diffusion_steps // val: int
    protenix_dtype         // val: string
    chai1_use_templates    // val: boolean
    binder_chains          // val: string
    target_chains          // val: string
    contact_cutoff         // val: float
    hbond_cutoff           // val: float
    clash_cutoff           // val: float

    main:
    ch_versions = Channel.empty()

    LIGANDMPNN(
        ch_structures,
        ch_mpnn_weights,
        mpnn_model_type,
        num_seq_per_pdb,
        mpnn_batch_size,
        mpnn_temperature,
        mpnn_seed,
        mpnn_output_mode,
        fixed_residues,
        chains_to_design
    )
    ch_versions = ch_versions.mix(LIGANDMPNN.out.versions)

    ch_pred_inputs = LIGANDMPNN.out.sequences_tsv.map { meta, seq_tsv ->
        [meta, [], seq_tsv]
    }

    if (predictor == 'alphafold3') {
        ALPHAFOLD3(
            ch_pred_inputs,
            ch_af3_model_dir,
            ch_af3_db_dir,
            target_sequence,
            ligand_ccd,
            ligand_smiles,
            af3_model_seeds,
            af3_run_data_pipeline,
            diffusion_samples
        )
        ch_pred_structures = ALPHAFOLD3.out.structures
        ch_pred_metrics_tsv = ALPHAFOLD3.out.metrics_tsv
        ch_versions = ch_versions.mix(ALPHAFOLD3.out.versions)
    } else if (predictor == 'protenix') {
        PROTENIX(
            ch_pred_inputs,
            ch_protenix_cache,
            protenix_model_name,
            target_sequence,
            ligand_ccd,
            ligand_smiles,
            af3_model_seeds,
            recycling_steps,
            protenix_diffusion_steps,
            diffusion_samples,
            use_msa_server,
            protenix_dtype
        )
        ch_pred_structures = PROTENIX.out.structures
        ch_pred_metrics_tsv = PROTENIX.out.metrics_tsv
        ch_versions = ch_versions.mix(PROTENIX.out.versions)
    } else if (predictor == 'chai1') {
        CHAI1(
            ch_pred_inputs,
            ch_chai1_cache,
            target_sequence,
            ligand_smiles,
            use_msa_server,
            chai1_use_templates,
            recycling_steps,
            diffusion_samples,
            mpnn_seed
        )
        ch_pred_structures = CHAI1.out.structures
        ch_pred_metrics_tsv = CHAI1.out.metrics_tsv
        ch_versions = ch_versions.mix(CHAI1.out.versions)
    } else if (predictor == 'boltz2') {
        BOLTZ2(
            ch_pred_inputs,
            ch_boltz_cache,
            target_sequence,
            ligand_ccd,
            ligand_smiles,
            predict_affinity,
            use_msa_server,
            use_potentials,
            recycling_steps,
            diffusion_samples,
            boltz_output_format
        )
        ch_pred_structures = BOLTZ2.out.structures
        ch_pred_metrics_tsv = BOLTZ2.out.metrics_tsv
        ch_versions = ch_versions.mix(BOLTZ2.out.versions)
    } else {
        error "Unsupported predictor '${predictor}' in SEQ_FOLD_EVAL. Expected 'boltz2', 'alphafold3', 'protenix', or 'chai1'."
    }

    PROTEIN_METRICS(
        ch_pred_structures,
        ch_reference,
        binder_chains,
        target_chains,
        contact_cutoff,
        hbond_cutoff,
        clash_cutoff
    )
    ch_versions = ch_versions.mix(PROTEIN_METRICS.out.versions)

    ch_all_tables = LIGANDMPNN.out.sequences_tsv
        .join(ch_pred_metrics_tsv, by: 0)
        .join(PROTEIN_METRICS.out.metrics_tsv, by: 0)
        .map { meta, mpnn_tsv, pred_tsv, struct_tsv ->
            [meta, [mpnn_tsv, pred_tsv, struct_tsv]]
        }

    emit:
    sequences_tsv         = LIGANDMPNN.out.sequences_tsv     // channel: [ val(meta), path(tsv) ]
    sequences_fasta       = LIGANDMPNN.out.sequences_fasta   // channel: [ val(meta), path(fasta) ]
    predicted_structures  = ch_pred_structures               // channel: [ val(meta), [ path(cif|pdb) ] ]
    predictor_metrics_tsv = ch_pred_metrics_tsv              // channel: [ val(meta), path(tsv) ]
    structure_metrics_tsv = PROTEIN_METRICS.out.metrics_tsv  // channel: [ val(meta), path(tsv) ]
    all_metric_tables     = ch_all_tables                    // channel: [ val(meta), [ mpnn_tsv, pred_tsv, struct_tsv ] ]
    versions              = ch_versions                      // channel: [ path(versions.yml) ]
}
