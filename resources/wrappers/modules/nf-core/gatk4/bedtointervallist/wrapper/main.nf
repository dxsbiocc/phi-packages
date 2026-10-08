#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { GATK4_BEDTOINTERVALLIST } from '../main.nf'


params.bed        = null
params.dict       = null
params.outdir     = null

workflow {
    bed_ch  = Channel.value([[id: 'test'], [file(params.bed, checkIfExists: true)]])
    dict_ch = Channel.value([[id: 'dict'], [file(params.dict, checkIfExists: true)]])

    GATK4_BEDTOINTERVALLIST(bed_ch, dict_ch)
}
