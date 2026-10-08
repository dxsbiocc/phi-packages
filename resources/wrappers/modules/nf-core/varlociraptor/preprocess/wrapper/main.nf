#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { VARLOCIRAPTOR_PREPROCESS } from '../main.nf'
include { VARLOCIRAPTOR_ESTIMATEALIGNMENTPROPERTIES } from '../../estimatealignmentproperties/main.nf'

params.bam        = null
params.bai        = null
params.fasta      = null
params.fai        = null
params.candidates = null
params.outdir     = null

workflow {
    est_in = Channel.value([[id: 'test_normal'], file(params.bam, checkIfExists: true), file(params.bai, checkIfExists: true), file(params.fasta, checkIfExists: true), file(params.fai, checkIfExists: true)])
    VARLOCIRAPTOR_ESTIMATEALIGNMENTPROPERTIES(est_in)
    pre_in = Channel.of([[id: 'test_normal'], file(params.bam, checkIfExists: true), file(params.bai, checkIfExists: true), file(params.candidates, checkIfExists: true),
                         file(params.fasta, checkIfExists: true), file(params.fai, checkIfExists: true)])
        .collect()
        .join(VARLOCIRAPTOR_ESTIMATEALIGNMENTPROPERTIES.out.alignment_properties_json)
        .map { meta, bam, bai, candidates, fasta, fai, alignment_json -> [meta, bam, bai, candidates, alignment_json, fasta, fai] }
    VARLOCIRAPTOR_PREPROCESS(pre_in)
}
