process FIGURES {
    cpus 2
    memory '8 GB'
    time '1h'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path qc_dir
    path ann_dir
    path state_dir
    path spatial_dir
    path posthoc_dir

    output:
    path "06_figures"

    script:
    """
    python ${projectDir}/notebooks/06_figures.py \\
        --results . \\
        --labels ${projectDir}/config/cluster_labels.yaml \\
        --outdir 06_figures
    """
}