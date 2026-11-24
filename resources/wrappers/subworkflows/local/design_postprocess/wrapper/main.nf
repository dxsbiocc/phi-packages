#!/usr/bin/env nextflow
// Agent-facing wrapper adapter for local/subworkflows/design-postprocess.
// Follows docs/contracts/wrapper.md.
nextflow.enable.dsl = 2

include { DESIGN_POSTPROCESS } from '../main.nf'

params.structures = null
params.extra_tables = null
params.reference = null
params.binder_chains = null
params.target_chains = null
params.contact_cutoff = null
params.hbond_cutoff = null
params.clash_cutoff = null
params.run_foldseek_cluster = null
params.foldseek_coverage = null
params.foldseek_tmscore_threshold = null
params.id_column = null
params.join_mode = null
params.expr = null
params.rank_by = null
params.rank_direction = null
params.group_by = null
params.max_per_group = null
params.top_n = null
params.outdir = null

workflow {
    def meta = [id: 'design_postprocess']
    def struct_files = file(params.structures, checkIfExists: true)
    def struct_list = struct_files instanceof List ? struct_files : [struct_files]

    def extra_files = params.extra_tables ? file(params.extra_tables, checkIfExists: true) : []
    def extra_list = extra_files instanceof List ? extra_files : (extra_files ? [extra_files] : [])
    def ref_file = params.reference ? file(params.reference, checkIfExists: true) : []

    DESIGN_POSTPROCESS(
        Channel.value([meta, struct_list]),
        Channel.value([meta, extra_list]),
        Channel.value(ref_file),
        params.binder_chains ?: '',
        params.target_chains ?: '',
        params.contact_cutoff,
        params.hbond_cutoff,
        params.clash_cutoff,
        params.run_foldseek_cluster,
        params.foldseek_coverage,
        params.foldseek_tmscore_threshold,
        params.id_column,
        params.join_mode,
        params.expr ?: '',
        params.rank_by ?: '',
        params.rank_direction,
        params.group_by ?: '',
        params.max_per_group,
        params.top_n
    )
}
