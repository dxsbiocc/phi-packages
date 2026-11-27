//
// Alternative Splicing Quantification Workflow (rMATS prep + SUPPA2)
//
// Composes vendored nf-core modules:
//   - RMATS_PREP (BAM -> *.rmats prep counts + read outcomes)
//   - SUPPA_GENERATEEVENTS (GTF -> *.ioe events)
//   - SUPPA_PSIPEREVENT (TPM expression + .ioe -> *.psi)
//   - SUPPA_DIFFSPLICE (condition1 vs condition2 -> *.dpsi + *.psivec)
//
// Outputs (*.rmats, *.dpsi) feed directly into bulk-rna-seq's rna_splicing tool.
//

include { RMATS_PREP           } from '../../../modules/nf-core/rmats/prep/main.nf'
include { SUPPA_GENERATEEVENTS } from '../../../modules/nf-core/suppa/generateevents/main.nf'
include { SUPPA_PSIPEREVENT    } from '../../../modules/nf-core/suppa/psiperevent/main.nf'
include { SUPPA_DIFFSPLICE     } from '../../../modules/nf-core/suppa/diffsplice/main.nf'

workflow RNASPLICE {
    take:
    ch_gtf             // channel: [ val(meta), path(gtf) ]
    ch_bams            // channel: [ val(meta), path(bam) ] (optional for rmats)
    ch_tx_cond1        // channel: [ val(meta), path(tpm_cond1) ] (optional for suppa)
    ch_tx_cond2        // channel: [ val(meta), path(tpm_cond2) ] (optional for suppa)
    run_rmats          // val: boolean
    run_suppa          // val: boolean
    read_length        // val: int
    suppa_event_type   // val: string (e.g. 'SE SS MX RI FL')
    suppa_method       // val: string ('empirical' | 'classical')
    condition1_name    // val: string
    condition2_name    // val: string

    main:
    ch_rmats = Channel.empty()
    ch_dpsi  = Channel.empty()
    ch_psi   = Channel.empty()

    if (run_rmats) {
        RMATS_PREP(
            ch_bams,
            ch_gtf.first(),
            read_length
        )
        ch_rmats = RMATS_PREP.out.rmats
    }

    if (run_suppa) {
        SUPPA_GENERATEEVENTS(
            ch_gtf,
            'ioe',
            true,
            suppa_event_type,
            'S',
            10,
            100
        )

        ch_events = SUPPA_GENERATEEVENTS.out.events.first()
        ch_expr_both = ch_tx_cond1
            .map { meta, tpm -> tuple(meta + [cond: condition1_name], tpm) }
            .mix(
                ch_tx_cond2.map { meta, tpm -> tuple(meta + [cond: condition2_name], tpm) }
            )

        SUPPA_PSIPEREVENT(
            ch_expr_both,
            ch_events,
            0
        )
        ch_psi = SUPPA_PSIPEREVENT.out.psi

        ch_cond1_pair = ch_tx_cond1.join(
            ch_psi.filter { meta, psi -> meta.cond == condition1_name }
                  .map { meta, psi -> tuple([id: meta.id], psi) }
        ).map { meta, tpm, psi -> tuple(meta, condition1_name, tpm, psi) }

        ch_cond2_pair = ch_tx_cond2.join(
            ch_psi.filter { meta, psi -> meta.cond == condition2_name }
                  .map { meta, psi -> tuple([id: meta.id], psi) }
        ).map { meta, tpm, psi -> tuple(meta, condition2_name, tpm, psi) }

        SUPPA_DIFFSPLICE(
            ch_events,
            ch_cond1_pair,
            ch_cond2_pair,
            suppa_method,
            20,
            0.05,
            false,
            true,
            0.05,
            true,
            false,
            false,
            0,
            0
        )
        ch_dpsi = SUPPA_DIFFSPLICE.out.dpsi
    }

    emit:
    rmats = ch_rmats
    psi   = ch_psi
    dpsi  = ch_dpsi
}
