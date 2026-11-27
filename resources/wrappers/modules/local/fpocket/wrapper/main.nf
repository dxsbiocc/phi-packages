nextflow.enable.dsl = 2

include { FPOCKET } from '../main.nf'

workflow {
    ch_input = Channel.of(
        tuple(
            [id: params.id ?: 'fpocket_scan'],
            file(params.receptor, checkIfExists: true)
        )
    )

    FPOCKET(
        ch_input,
        params.min_alpha_spheres ?: 35,
        params.box_padding ?: 6.0
    )
}
