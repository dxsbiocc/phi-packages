#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { SNPSIFT_DBNSFP } from '../main.nf'


params.vcf        = null
params.tbi        = null
params.database   = null
params.database_tbi = null
params.outdir     = null

workflow {
    vcf_ch = Channel.value([[id: 'sample', single_end: false], file(params.vcf, checkIfExists: true), file(params.tbi, checkIfExists: true)])
    db_ch  = Channel.value([[id: 'databases'], file(params.database, checkIfExists: true), file(params.database_tbi, checkIfExists: true)])

    SNPSIFT_DBNSFP(vcf_ch, db_ch)
}
