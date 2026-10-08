#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { ELPREP_FILTER } from '../main.nf'


params.bam        = null
params.bai        = null
params.target_regions = null
params.elfasta    = null
params.elsites    = null
params.outdir     = null

workflow {
    in_ch = Channel.value([[id: 'sample', single_end: false], file(params.bam, checkIfExists: true), file(params.bai, checkIfExists: true),
                           file(params.target_regions, checkIfExists: true), [], [], []])
    elfasta_ch = Channel.value([[id: 'elfasta'], file(params.elfasta, checkIfExists: true)])
    elsites_ch = Channel.value([[id: 'sites'], file(params.elsites, checkIfExists: true)])

    ELPREP_FILTER(in_ch, [[], []], elfasta_ch, elsites_ch, true, true, false, false, false)
}
