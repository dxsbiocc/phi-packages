nextflow.enable.dsl = 2

include { GNINA } from '../main.nf'

workflow {
    ch_autobox = params.autobox_ligand ? file(params.autobox_ligand, checkIfExists: true) : []
    ch_input = Channel.of(
        tuple(
            [id: params.id ?: 'gnina_screen'],
            file(params.receptor, checkIfExists: true),
            file(params.ligands, checkIfExists: true),
            ch_autobox
        )
    )

    GNINA(
        ch_input,
        params.center_x ?: 0.0,
        params.center_y ?: 0.0,
        params.center_z ?: 0.0,
        params.size_x ?: 22.0,
        params.size_y ?: 22.0,
        params.size_z ?: 22.0,
        params.cnn_scoring ?: 'rescore',
        params.exhaustiveness ?: 8,
        params.num_modes ?: 9,
        params.cpu_threads ?: 4
    )
}
