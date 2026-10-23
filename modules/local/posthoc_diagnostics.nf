process POSTHOC_DIAGNOSTICS {
    cpus 2
    memory '8 GB'
    time '2h'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path qc_dir

    output:
    path "05_posthoc"

    script:
    """
    python ${projectDir}/notebooks/05_posthoc_diagnostics.py \\
        --h5ad_dir ${qc_dir}/h5ad \\
        --outdir 05_posthoc
    """
}