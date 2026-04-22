process SUBCLONAL_INTERPRETATION {
    tag "$meta.id"
    label "process_low"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5a/5abea2c371bb95e1db0081b89849fcd89015844fd2deb09fda56aefeb87a45f7/data':
        'community.wave.seqera.io/library/r-ggplot2_r-ggrepel_r-tidyverse:29f625041d233719' }"

    input:
    tuple val(meta), path(mutation_tables), path(results_sigprofiler)
    // one element in the tuple for dataset
    // { (meta1, [mutation tables list], [sigprofiler folders list]),
    //   (meta2, [mutation tables list], [sigprofiler folders list]), ... }
    // mutation_tables -> list of tables output from `prepare_cluster`
    // results_sigprofiler -> path to sigprofiler folder output from `assign_cluster`

    output:
    tuple val(meta), path("*.pdf"), emit: report_pdf

    when:
    task.ext.when == null || task.ext.when

    script:
    template "main_script.R"

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    """
}
