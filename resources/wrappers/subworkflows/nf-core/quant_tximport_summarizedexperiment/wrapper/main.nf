#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored subworkflow at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1. Includes
// the subworkflow only — never another wrapper (section 6).
nextflow.enable.dsl = 2

include { QUANT_TXIMPORT_SUMMARIZEDEXPERIMENT } from '../main.nf'

params.quants        = null
params.gtf           = null
params.samplesheet   = null
params.quant_type    = 'kallisto'
params.id            = 'gene_id'
params.extra         = 'gene_name'
params.skip_merge    = false
params.outdir        = null

workflow {
    quants_ch = Channel
        .fromPath(params.quants, checkIfExists: true, type: 'dir')
        .map { dir -> [[id: dir.name], dir] }

    gtf_ch = Channel.of(file(params.gtf, checkIfExists: true))

    samplesheet_ch = Channel.value([[id: 'samplesheet'], file(params.samplesheet, checkIfExists: true)])

    QUANT_TXIMPORT_SUMMARIZEDEXPERIMENT(
        samplesheet_ch,
        quants_ch,
        gtf_ch,
        params.id,
        params.extra,
        params.quant_type,
        params.skip_merge
    )
}
