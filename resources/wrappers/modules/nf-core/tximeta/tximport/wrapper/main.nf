#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
// Composes custom/tx2gene directly (vendored module, not another wrapper)
// to build the transcript-to-gene map tximport needs.
nextflow.enable.dsl = 2

include { CUSTOM_TX2GENE } from '../../../custom/tx2gene/main.nf'
include { TXIMETA_TXIMPORT } from '../main.nf'

params.gtf        = null
params.quants     = null
params.quant_type = 'kallisto'
params.id         = 'gene_id'
params.extra      = 'gene_name'
params.outdir     = null

workflow {
    gtf_ch = Channel
        .fromPath(params.gtf, checkIfExists: true)
        .map { gtf -> [[id: 'test'], gtf] }

    quants_ch = Channel
        .fromPath(params.quants, checkIfExists: true, type: 'dir')
        .collect()
        .map { dirs -> [[id: 'test'], dirs] }

    CUSTOM_TX2GENE(gtf_ch, quants_ch, params.quant_type, params.id, params.extra)

    TXIMETA_TXIMPORT(quants_ch, CUSTOM_TX2GENE.out.tx2gene, params.quant_type)
}
