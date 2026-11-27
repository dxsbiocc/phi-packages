nextflow.enable.dsl = 2

include { AUTODOCK_VINA } from '../main.nf'

workflow {
    ch_input = Channel.of(
        tuple(
            [id: params.id ?: 'vina_screen'],
            file(params.receptor, checkIfExists: true),
            file(params.ligands_dir, checkIfExists: true)
        )
    )

    AUTODOCK_VINA(
        ch_input,
        params.center_x,
        params.center_y,
        params.center_z,
        params.size_x ?: 22.0,
        params.size_y ?: 22.0,
        params.size_z ?: 22.0,
        params.exhaustiveness ?: 8,
        params.num_modes ?: 9,
        params.cpu_per_ligand ?: 4
    )
}
