#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { MULTIQC } from '../main.nf'

params.logs   = null
params.outdir = null

workflow {
    logs_ch = Channel
        .fromPath(params.logs, checkIfExists: true)
        .collect()
        .map { files -> [[id: 'multiqc'], files, [], [], [], []] }

    MULTIQC(logs_ch)
}
