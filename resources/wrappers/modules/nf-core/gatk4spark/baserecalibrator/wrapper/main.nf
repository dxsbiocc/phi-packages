#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { GATK4SPARK_BASERECALIBRATOR } from '../main.nf'


params.bam        = null
params.bai        = null
params.fasta      = null
params.fai        = null
params.dict       = null
params.known_sites = null
params.known_sites_tbi = null
params.outdir     = null

workflow {
    in_ch = Channel.of([[id: 'sample'], file(params.bam, checkIfExists: true), file(params.bai, checkIfExists: true), []])

    GATK4SPARK_BASERECALIBRATOR(in_ch, Channel.of(file(params.fasta, checkIfExists: true)), Channel.of(file(params.fai, checkIfExists: true)),
                                Channel.of(file(params.dict, checkIfExists: true)), Channel.of(file(params.known_sites, checkIfExists: true)), Channel.of(file(params.known_sites_tbi, checkIfExists: true)))
}
