#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { SYLPHTAX_MERGE } from '../main.nf'
include { SYLPH_PROFILE } from '../../../sylph/profile/main.nf'
include { SYLPHTAX_TAXPROF } from '../../taxprof/main.nf'

params.reads1     = null
params.reads2     = null
params.database   = null
params.taxonomy   = null
params.data_type  = 'relative_abundance'
params.outdir     = null

workflow {
    reads_ch = Channel.value([[id: 'sample', single_end: false], [file(params.reads1, checkIfExists: true), file(params.reads2, checkIfExists: true)]])

    SYLPH_PROFILE(reads_ch, file(params.database, checkIfExists: true))
    SYLPHTAX_TAXPROF(SYLPH_PROFILE.out.profile_out, file(params.taxonomy, checkIfExists: true))
    SYLPHTAX_MERGE(SYLPHTAX_TAXPROF.out.taxprof_output, params.data_type)
}
