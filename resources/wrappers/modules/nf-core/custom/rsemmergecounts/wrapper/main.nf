#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { CUSTOM_RSEMMERGECOUNTS } from '../main.nf'

params.genes_results    = null
params.isoforms_results = null
params.outdir           = null

workflow {
    genes_ch = Channel
        .fromPath(params.genes_results, checkIfExists: true)
        .collect()
        .map { files -> [[id: 'rsem_merged'], files] }

    isoforms_ch = Channel
        .fromPath(params.isoforms_results, checkIfExists: true)
        .collect()

    CUSTOM_RSEMMERGECOUNTS(genes_ch, isoforms_ch)
}
