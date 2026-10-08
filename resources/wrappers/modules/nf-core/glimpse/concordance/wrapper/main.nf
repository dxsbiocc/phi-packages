#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { GLIMPSE_CONCORDANCE } from '../main.nf'


params.estimate   = null
params.estimate_csi = null
params.freq       = null
params.freq_csi   = null
params.truth      = null
params.truth_csi  = null
params.region     = 'chr21'
params.outdir     = null

workflow {
    in_ch = Channel.value([[id: 'NA12878'], file(params.estimate, checkIfExists: true), file(params.estimate_csi, checkIfExists: true),
                           file(params.freq, checkIfExists: true), file(params.freq_csi, checkIfExists: true),
                           file(params.truth, checkIfExists: true), file(params.truth_csi, checkIfExists: true), params.region])

    GLIMPSE_CONCORDANCE(in_ch, 0.7, 3, [])
}
