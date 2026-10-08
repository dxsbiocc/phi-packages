#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { PLINK2_PMERGE } from '../main.nf'
include { PLINK2_FILTER as PLINK2_FILTER_PART1 } from '../../filter/main.nf'
include { PLINK2_FILTER as PLINK2_FILTER_PART2 } from '../../filter/main.nf'

params.pgen       = null
params.pvar       = null
params.psam       = null
params.outdir     = null

workflow {
    part1_in = Channel.value([[id: 'part1'], file(params.pgen, checkIfExists: true), file(params.pvar, checkIfExists: true), file(params.psam, checkIfExists: true)])
    part2_in = Channel.value([[id: 'part2'], file(params.pgen, checkIfExists: true), file(params.pvar, checkIfExists: true), file(params.psam, checkIfExists: true)])
    PLINK2_FILTER_PART1(part1_in)
    PLINK2_FILTER_PART2(part2_in)

    merge_a = PLINK2_FILTER_PART1.out.pgen.join(PLINK2_FILTER_PART1.out.pvar).join(PLINK2_FILTER_PART1.out.psam)
        .map { meta, pgen, pvar, psam -> [[id: 'test'], pgen, pvar, psam] }
    merge_b = PLINK2_FILTER_PART2.out.pgen.join(PLINK2_FILTER_PART2.out.pvar).join(PLINK2_FILTER_PART2.out.psam)
        .map { meta, pgen, pvar, psam -> [[id: 'test'], pgen, pvar, psam] }

    PLINK2_PMERGE(merge_a, merge_b)

}
