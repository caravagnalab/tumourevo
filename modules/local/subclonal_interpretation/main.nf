process SUBCLONAL_INTERPRETATION {
    tag "$meta.id"
    label "process_low"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/53/53da5a4364ce10f58c67142c1b4203f19eeff99b533c3e3a310d0a3b37fb46e8/data':
        'community.wave.seqera.io/library/r-ggforce_r-ggnewscale_r-ggraph_r-ggrepel_pruned:5408b55f63d978d5' }"

    input:
    tuple val(meta), path(mutation_tables), path(results_sigprofiler), path(fits, stageAs: "fits/*"), path(joint_table, stageAs: "joint/*")

    output:
    tuple val(meta), path("*_scores_clusters.pdf"), emit: report_score
    tuple val(meta), path("*_signature_clusters.pdf"), emit: report_signature
    tuple val(meta), path("*_scores.rds"), emit: rds_score
    tuple val(meta), path("*_table_signature.rds"), emit: rds_signature
    tuple val(meta), path("*_scores.tsv"), emit: tsv_score
    tuple val(meta), path("*_sample_quality.tsv"), emit: tsv_sample_quality

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    def tools = args!="" && args.tools ? "$args.tools" : ""

    template "main_script.R"

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}_scores_clusters.pdf
    touch ${prefix}_signature_clusters.pdf
    touch ${prefix}_scores.rds
    touch ${prefix}_table_signature.rds
    touch ${prefix}_scores.tsv
    touch ${prefix}_sample_quality.tsv
    """
}
