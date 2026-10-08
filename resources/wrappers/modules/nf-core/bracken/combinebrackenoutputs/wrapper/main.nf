#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { BRACKEN_COMBINEBRACKENOUTPUTS } from '../main.nf'
include { UNTAR } from '../../../untar/main.nf'
include { KRAKEN2_KRAKEN2 } from '../../../kraken2/kraken2/main.nf'
include { BRACKEN_BRACKEN } from '../../bracken/main.nf'

params.reads1     = null
params.reads2     = null
params.db         = null
params.outdir     = null

workflow {
    UNTAR(Channel.value([[id: 'bracken_db'], file(params.db, checkIfExists: true)]))
    db_dir_ch = UNTAR.out.untar.map { meta, dir -> dir }
    reads_ch = Channel.of([[id: 'sample1', single_end: true], file(params.reads1, checkIfExists: true)],
                          [[id: 'sample2', single_end: true], file(params.reads2, checkIfExists: true)])

    KRAKEN2_KRAKEN2(reads_ch, db_dir_ch, false, false)
    BRACKEN_BRACKEN(KRAKEN2_KRAKEN2.out.report, db_dir_ch)
    BRACKEN_COMBINEBRACKENOUTPUTS(BRACKEN_BRACKEN.out.reports.map { it[1] }.collect().map { [[id: 'db'], it] })
}
