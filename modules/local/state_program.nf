process STATE_PROGRAM {
    cpus 2
    memory '8 GB'
    time '1h'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path ann_dir

    output:
    path "03_state"

    script:
    """
    python ${projectDir}/notebooks/03_state_program.py \\
        --h5ad_dir ${ann_dir} \\
        --labels ${projectDir}/config/cluster_labels.yaml \\
        --outdir 03_state
    """
}