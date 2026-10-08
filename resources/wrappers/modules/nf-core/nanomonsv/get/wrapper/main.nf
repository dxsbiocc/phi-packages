#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { NANOMONSV_GET } from '../main.nf'
include { NANOMONSV_PARSE as NANOMONSV_PARSE_TUMOR } from '../../parse/main.nf'
include { NANOMONSV_PARSE as NANOMONSV_PARSE_CONTROL } from '../../parse/main.nf'

params.tumor_bam  = null
params.tumor_bai  = null
params.control_bam = null
params.control_bai = null
params.fasta      = null
params.fai        = null
params.outdir     = null

workflow {
    tumor_bam_ch = Channel.value([[id: 'tumor'], file(params.tumor_bam, checkIfExists: true), file(params.tumor_bai, checkIfExists: true)])
    control_bam_ch = Channel.value([[id: 'control'], file(params.control_bam, checkIfExists: true), file(params.control_bai, checkIfExists: true)])
    NANOMONSV_PARSE_TUMOR(tumor_bam_ch)
    NANOMONSV_PARSE_CONTROL(control_bam_ch)

    tumor_parse = NANOMONSV_PARSE_TUMOR.out.insertions.combine(NANOMONSV_PARSE_TUMOR.out.insertions_index, by: 0)
        .combine(NANOMONSV_PARSE_TUMOR.out.deletions, by: 0).combine(NANOMONSV_PARSE_TUMOR.out.deletions_index, by: 0)
        .combine(NANOMONSV_PARSE_TUMOR.out.rearrangements, by: 0).combine(NANOMONSV_PARSE_TUMOR.out.rearrangements_index, by: 0)
        .combine(NANOMONSV_PARSE_TUMOR.out.bp_info, by: 0).combine(NANOMONSV_PARSE_TUMOR.out.bp_info_index, by: 0)
        .map { meta, ins, ins_tbi, del, del_tbi, rea, rea_tbi, bp, bp_tbi -> [meta, [ins, ins_tbi, del, del_tbi, rea, rea_tbi, bp, bp_tbi]] }
    control_parse = NANOMONSV_PARSE_CONTROL.out.insertions.combine(NANOMONSV_PARSE_CONTROL.out.insertions_index, by: 0)
        .combine(NANOMONSV_PARSE_CONTROL.out.deletions, by: 0).combine(NANOMONSV_PARSE_CONTROL.out.deletions_index, by: 0)
        .combine(NANOMONSV_PARSE_CONTROL.out.rearrangements, by: 0).combine(NANOMONSV_PARSE_CONTROL.out.rearrangements_index, by: 0)
        .combine(NANOMONSV_PARSE_CONTROL.out.bp_info, by: 0).combine(NANOMONSV_PARSE_CONTROL.out.bp_info_index, by: 0)
        .map { meta, ins, ins_tbi, del, del_tbi, rea, rea_tbi, bp, bp_tbi -> [meta, [ins, ins_tbi, del, del_tbi, rea, rea_tbi, bp, bp_tbi]] }

    in0 = tumor_bam_ch.combine(tumor_parse, by: 0).map { meta, bam, bai, parse_files -> [meta, bam, bai, parse_files] }
    in1 = control_bam_ch.combine(control_parse, by: 0).map { meta, bam, bai, parse_files -> [meta, bam, bai, parse_files] }
    fasta_ch = Channel.value([[id: 'genome'], file(params.fasta, checkIfExists: true), file(params.fai, checkIfExists: true)])

    NANOMONSV_GET(in0, in1, fasta_ch, [], [], [])
}
