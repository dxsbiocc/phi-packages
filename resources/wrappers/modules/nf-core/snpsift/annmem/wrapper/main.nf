#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { SNPSIFT_ANNMEM } from '../main.nf'
include { SNPSIFT_ANNMEMCREATE } from '../../annmemcreate/main.nf'

params.vcf        = null
params.tbi        = null
params.db_vcf     = null
params.db_tbi     = null
params.fields     = 'DP,VDB'
params.outdir     = null

workflow {
    vcf_ch = Channel.value([[id: 'sample'], file(params.vcf, checkIfExists: true), file(params.tbi, checkIfExists: true)])
    SNPSIFT_ANNMEMCREATE(Channel.value([[id: 'db'], file(params.vcf, checkIfExists: true), file(params.tbi, checkIfExists: true), params.fields.tokenize(',')]))
    db_ch = SNPSIFT_ANNMEMCREATE.out.database.map { meta, vardb -> [file(params.db_vcf, checkIfExists: true), file(params.db_tbi, checkIfExists: true), vardb, [], []] }

    SNPSIFT_ANNMEM(vcf_ch, db_ch)
}
