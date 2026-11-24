//
// Protein Design Post-processing Subworkflow
//
// Evaluates designed or predicted protein structures with PROTEIN_METRICS
// (Kabsch RMSD, interface contacts, putative H-bonds, ligand contacts, clashes,
// radius of gyration), optionally runs FOLDSEEK (easy-cluster) to assign
// structural `cluster_id` groups, then merges all upstream, structural, and
// cluster tables on design_id and runs Chaperone-compatible boolean DSL
// filtering, ranking, and diversity selection via DESIGN_FILTER.
//

include { PROTEIN_METRICS } from '../../../modules/local/protein_metrics/main.nf'
include { FOLDSEEK        } from '../../../modules/local/foldseek/main.nf'
include { DESIGN_FILTER   } from '../../../modules/local/design_filter/main.nf'

workflow DESIGN_POSTPROCESS {
    take:
    ch_structures              // channel: [ val(meta), [ path(structure), ... ] ]
    ch_extra_tables            // channel: [ val(meta), [ path(table), ... ] ]  (may be empty list per meta)
    ch_reference               // channel: path(reference) or []
    binder_chains              // val: string
    target_chains              // val: string
    contact_cutoff             // val: float
    hbond_cutoff               // val: float
    clash_cutoff               // val: float
    run_foldseek_cluster       // val: boolean
    foldseek_coverage          // val: float
    foldseek_tmscore_threshold // val: float
    id_column                  // val: string
    join_mode                  // val: string
    expr                       // val: string
    rank_by                    // val: string
    rank_direction             // val: string
    group_by                   // val: string
    max_per_group              // val: int
    top_n                      // val: int

    main:
    ch_versions = Channel.empty()

    PROTEIN_METRICS(
        ch_structures,
        ch_reference,
        binder_chains,
        target_chains,
        contact_cutoff,
        hbond_cutoff,
        clash_cutoff
    )
    ch_versions = ch_versions.mix(PROTEIN_METRICS.out.versions)

    if (run_foldseek_cluster) {
        FOLDSEEK(
            ch_structures,
            [],
            'cluster',
            0.01,
            foldseek_coverage,
            foldseek_tmscore_threshold,
            0.5
        )
        ch_foldseek_tsv = FOLDSEEK.out.metrics_tsv
        ch_versions = ch_versions.mix(FOLDSEEK.out.versions)

        ch_filter_inputs = PROTEIN_METRICS.out.metrics_tsv
            .join(ch_foldseek_tsv, by: 0)
            .join(ch_extra_tables, by: 0, remainder: true)
            .map { meta, metrics_tsv, foldseek_tsv, extra_tables ->
                def all_tables = []
                if (extra_tables) {
                    all_tables.addAll(extra_tables instanceof List ? extra_tables : [extra_tables])
                }
                if (metrics_tsv) {
                    all_tables.add(metrics_tsv)
                }
                if (foldseek_tsv) {
                    all_tables.add(foldseek_tsv)
                }
                [meta, all_tables]
            }
    } else {
        ch_foldseek_tsv = Channel.empty()

        ch_filter_inputs = PROTEIN_METRICS.out.metrics_tsv
            .join(ch_extra_tables, by: 0, remainder: true)
            .map { meta, metrics_tsv, extra_tables ->
                def all_tables = []
                if (extra_tables) {
                    all_tables.addAll(extra_tables instanceof List ? extra_tables : [extra_tables])
                }
                if (metrics_tsv) {
                    all_tables.add(metrics_tsv)
                }
                [meta, all_tables]
            }
    }

    def effective_group_by = group_by ?: (run_foldseek_cluster ? 'cluster_id' : '')
    def effective_max_per_group = (max_per_group && (max_per_group as int) > 0)
        ? (max_per_group as int)
        : (run_foldseek_cluster ? 1 : 0)

    DESIGN_FILTER(
        ch_filter_inputs,
        id_column,
        join_mode,
        expr,
        rank_by,
        rank_direction,
        effective_group_by,
        effective_max_per_group,
        top_n
    )
    ch_versions = ch_versions.mix(DESIGN_FILTER.out.versions)

    emit:
    metrics_tsv  = PROTEIN_METRICS.out.metrics_tsv   // channel: [ val(meta), path(tsv) ]
    metrics_json = PROTEIN_METRICS.out.metrics_json  // channel: [ val(meta), path(json) ]
    foldseek_tsv = ch_foldseek_tsv                   // channel: [ val(meta), path(tsv) ] (empty when disabled)
    merged_tsv   = DESIGN_FILTER.out.merged_tsv      // channel: [ val(meta), path(tsv) ]
    filtered_tsv = DESIGN_FILTER.out.filtered_tsv    // channel: [ val(meta), path(tsv) ]
    summary_json = DESIGN_FILTER.out.summary_json    // channel: [ val(meta), path(json) ]
    versions     = ch_versions                       // channel: [ path(versions.yml) ]
}
