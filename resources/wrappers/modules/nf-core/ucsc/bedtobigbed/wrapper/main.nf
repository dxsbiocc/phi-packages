#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { UCSC_BEDTOBIGBED } from '../main.nf'

params.bed    = null
params.sizes  = null
params.outdir = null

workflow {
    bed_ch   = Channel.value([[id: 'test'], [file(params.bed, checkIfExists: true)]])
    sizes_ch = Channel.value(file(params.sizes, checkIfExists: true))
    autosql_ch = Channel.value([])

    UCSC_BEDTOBIGBED(bed_ch, sizes_ch, autosql_ch)
}
