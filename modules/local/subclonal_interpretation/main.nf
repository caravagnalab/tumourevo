process SUBCLONAL_INTERPRETATION {
    tag "$meta.id"
    label "process_low"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5a/5abea2c371bb95e1db0081b89849fcd89015844fd2deb09fda56aefeb87a45f7/data':
        'community.wave.seqera.io/library/r-ggplot2_r-ggrepel_r-tidyverse:29f625041d233719' }"

    input:
    tuple val(meta), path(mutation_tables), path(results_sigprofiler)

    output:
    tuple val(meta), path("*.pdf"), emit: report_pdf
    tuple val(meta), path("*.rds"), emit: rds

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    template "main_script.R"

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}_scores_clusters.pdf
    touch ${prefix}_signature_clusters.pdf
    """
}
