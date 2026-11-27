nextflow.enable.dsl = 2

include { RNASPLICE } from '../main.nf'

workflow {
    def run_id = params.id ?: 'rnasplice_run'
    ch_gtf = Channel.of(tuple([id: run_id], file(params.gtf, checkIfExists: true)))

    ch_bams = params.bam_samplesheet
        ? Channel.fromPath(params.bam_samplesheet, checkIfExists: true)
            .splitCsv(header: true)
            .map { row -> tuple([id: row.sample], file(row.bam, checkIfExists: true)) }
        : Channel.empty()

    ch_tx_cond1 = params.tx_tpm_cond1
        ? Channel.of(tuple([id: run_id], file(params.tx_tpm_cond1, checkIfExists: true)))
        : Channel.empty()

    ch_tx_cond2 = params.tx_tpm_cond2
        ? Channel.of(tuple([id: run_id], file(params.tx_tpm_cond2, checkIfExists: true)))
        : Channel.empty()

    def do_rmats = (params.run_rmats != false) && (params.bam_samplesheet != null)
    def do_suppa = (params.run_suppa != false) && (params.tx_tpm_cond1 != null) && (params.tx_tpm_cond2 != null)

    RNASPLICE(
        ch_gtf,
        ch_bams,
        ch_tx_cond1,
        ch_tx_cond2,
        do_rmats,
        do_suppa,
        params.read_length ?: 150,
        params.suppa_event_type ?: 'SE SS MX RI FL',
        params.suppa_method ?: 'empirical',
        params.condition1_name ?: 'cond1',
        params.condition2_name ?: 'cond2'
    )
}
