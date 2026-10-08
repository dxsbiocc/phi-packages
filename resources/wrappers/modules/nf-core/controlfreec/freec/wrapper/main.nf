#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { CONTROLFREEC_FREEC } from '../main.nf'
include { UNTAR } from '../../../untar/main.nf'

params.normal     = null
params.tumor      = null
params.fasta      = null
params.fai        = null
params.known_snps = null
params.known_snps_tbi = null
params.chromosomes = null
params.target_bed = null
params.outdir     = null

workflow {
    UNTAR(Channel.value([[], file(params.chromosomes, checkIfExists: true)]))
    in_ch = Channel.value([[id: 'sample', single_end: false, sex: 'XX'], file(params.normal, checkIfExists: true), file(params.tumor, checkIfExists: true), [], [], [], []])
    CONTROLFREEC_FREEC(in_ch, file(params.fasta, checkIfExists: true), file(params.fai, checkIfExists: true), [],
                       file(params.known_snps, checkIfExists: true), file(params.known_snps_tbi, checkIfExists: true),
                       UNTAR.out.untar.map { it[1] }, [], file(params.target_bed, checkIfExists: true), [])
}
