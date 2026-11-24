#!/usr/bin/env nextflow
// Agent-facing wrapper adapter for local/modules/design-filter.
// Follows docs/contracts/wrapper.md.
nextflow.enable.dsl = 2

include { DESIGN_FILTER } from '../main.nf'

params.tables = null
params.id_column = null
params.join_mode = null
params.expr = null
params.rank_by = null
params.rank_direction = null
params.group_by = null
params.max_per_group = null
params.top_n = null
params.args = null
params.outdir = null

workflow {
    def table_files = file(params.tables, checkIfExists: true)
    def table_list = table_files instanceof List ? table_files : [table_files]

    DESIGN_FILTER(
        Channel.value([[id: 'design_filter'], table_list]),
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
