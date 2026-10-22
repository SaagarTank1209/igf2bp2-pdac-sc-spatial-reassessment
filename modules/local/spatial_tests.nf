process SPATIAL_TESTS {
    cpus 2
    memory '8 GB'
    time '2h'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    val h5ad_dir

    output:
    path "04_spatial"

    script:
    """
    python ${projectDir}/notebooks/04_spatial_tests.py \\
        --h5ad_dir ${h5ad_dir} \\
        --outdir 04_spatial
    """
}