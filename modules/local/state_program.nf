process STATE_PROGRAM {
    cpus 2
    memory '8 GB'
    time '1h'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    val h5ad_dir

    output:
    path "03_state"

    script:
    """
    python ${projectDir}/notebooks/03_state_program.py \\
        --h5ad_dir ${h5ad_dir} \\
        --labels ${projectDir}/config/cluster_labels.yaml \\
        --outdir 03_state
    """
}