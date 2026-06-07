process PLOT_DECONVOLUTION {
    tag "$meta.id"
    label "process_low"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/53/53da5a4364ce10f58c67142c1b4203f19eeff99b533c3e3a310d0a3b37fb46e8/data':
        'community.wave.seqera.io/library/r-ggforce_r-ggnewscale_r-ggraph_r-ggrepel_pruned:5408b55f63d978d5' }"

    input:
    tuple val(meta), path(rds_mobster), path(rds_viber), path(pyclone_table), path(pyclone_fit)

    output:
    tuple val(meta), path("*.pdf"), optional: true, emit: plot_pdf
    tuple val(meta), path("*.png"), optional: true, emit: plot_png

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    def tools = args!="" && args.tools ? "$args.tools" : ""

    template "main_script.R"

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.pdf
    """
}
