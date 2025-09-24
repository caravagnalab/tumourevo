process PYCLONEVI {
    tag "$meta.id"
    label "process_low"
    label "error_ignore"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pyclone-vi:0.1.6--pyhdfd78af_0' :
        'biocontainers/pyclone-vi:0.1.6--pyhdfd78af_0' }"

    input:
        tuple val(meta), path(rds_join), val(tumour_samples)

    output:
        tuple val(meta), path("*_cluster_table.csv"), emit: ctree_input
        tuple val(meta), path("*.tsv"), emit: pyclone_input
        tuple val(meta), path("*_all_fits.h5"), emit: pyclone_all_fits
        tuple val(meta), path("*_best_fit.txt"), emit: pyclone_best_fit
        path "versions.yml", emit: versions

    script:
    template "main_script.py"
}
