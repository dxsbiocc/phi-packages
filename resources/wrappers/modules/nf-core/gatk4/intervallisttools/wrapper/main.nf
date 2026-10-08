#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// Composes the vendored gatk4/bedtointervallist module (../../bedtointervallist/main.nf).
nextflow.enable.dsl = 2

include { GATK4_INTERVALLISTTOOLS } from '../main.nf'
include { GATK4_BEDTOINTERVALLIST } from '../../bedtointervallist/main.nf'

params.bed        = null
params.dict       = null
params.outdir     = null

workflow {
    bed_ch  = Channel.value([[id: 'test'], [file(params.bed, checkIfExists: true)]])
    dict_ch = Channel.value([[id: 'dict'], [file(params.dict, checkIfExists: true)]])

    GATK4_BEDTOINTERVALLIST(bed_ch, dict_ch)
    GATK4_INTERVALLISTTOOLS(GATK4_BEDTOINTERVALLIST.out.interval_list)
}
