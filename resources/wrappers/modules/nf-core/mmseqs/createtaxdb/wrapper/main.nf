#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { MMSEQS_CREATETAXDB } from '../main.nf'
include { UNTAR as UNTAR_MMSEQS_DB } from '../../../untar/main.nf'
include { UNTAR as UNTAR_TAXDUMP } from '../../../untar/main.nf'

params.db_archive = null
params.taxdump    = null
params.tax_mapping = null
params.outdir     = null

workflow {
    UNTAR_MMSEQS_DB(Channel.value([[id: 'test'], file(params.db_archive, checkIfExists: true)]))
    UNTAR_TAXDUMP(Channel.value([[id: 'taxdump'], file(params.taxdump, checkIfExists: true)]))

    MMSEQS_CREATETAXDB(UNTAR_MMSEQS_DB.out.untar, UNTAR_TAXDUMP.out.untar, Channel.value([[], file(params.tax_mapping, checkIfExists: true)]))
}
