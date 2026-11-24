//
// Protein Design Post-processing Subworkflow
//
// Evaluates designed or predicted protein structures with PROTEIN_METRICS
// (Kabsch RMSD, interface contacts, putative H-bonds, ligand contacts, clashes,
// radius of gyration), then merges all upstream and structural metric tables on
// design_id and runs Chaperone-compatible boolean DSL filtering, ranking, and
// per-group diversity selection via DESIGN_FILTER.
//

include { PROTEIN_METRICS } from '../../../modules/local/protein_metrics/main.nf'
include { DESIGN_FILTER   } from '../../../modules/local/design_filter/main.nf'

workflow DESIGN_POSTPROCESS {
    take:
    ch_structures     // channel: [ val(meta), [ path(structure), ... ] ]
    ch_extra_tables   // channel: [ val(meta), [ path(table), ... ] ]  (may be empty list per meta)
    ch_reference      // channel: path(reference) or []
    binder_chains     // val: string
    target_chains     // val: string
    contact_cutoff    // val: float
    hbond_cutoff      // val: float
    clash_cutoff      // val: float
    id_column         // val: string
    join_mode         // val: string
    expr              // val: string
    rank_by           // val: string
    rank_direction    // val: string
    group_by          // val: string
    max_per_group     // val: int
    top_n             // val: int

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

    DESIGN_FILTER(
        ch_filter_inputs,
        id_column,
        join_mode,
        expr,
        rank_by,
        rank_direction,
        group_by,
        max_per_group,
        top_n
    )
    ch_versions = ch_versions.mix(DESIGN_FILTER.out.versions)

    emit:
    metrics_tsv  = PROTEIN_METRICS.out.metrics_tsv   // channel: [ val(meta), path(tsv) ]
    metrics_json = PROTEIN_METRICS.out.metrics_json  // channel: [ val(meta), path(json) ]
    merged_tsv   = DESIGN_FILTER.out.merged_tsv      // channel: [ val(meta), path(tsv) ]
    filtered_tsv = DESIGN_FILTER.out.filtered_tsv    // channel: [ val(meta), path(tsv) ]
    summary_json = DESIGN_FILTER.out.summary_json    // channel: [ val(meta), path(json) ]
    versions     = ch_versions                       // channel: [ path(versions.yml) ]
}
