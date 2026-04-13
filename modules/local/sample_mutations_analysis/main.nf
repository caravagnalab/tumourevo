process SAMPLE_MUTATIONS_ANALYSIS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ca/ca9bc8ca16e42be96ad41d68b66f3c9645ed4e037c74ae9b4e236548d4599271/data':
        'community.wave.seqera.io/library/bioconductor-complexheatmap_r-dplyr_r-ggplot2_r-patchwork_r-tidyr:653ce8507801e2bb' }"
        
    input:
    tuple val(meta), path(snv_rds), val(tumour_sample)

    output:
    tuple val(meta), path("*_tmb.rds"), emit: tmb_rds 
    tuple val(meta), path("*_vaf_chr_plot.rds"), path('*_chr_mut.rds'), path('*_consequence_mut_plot.rds'), path('*_mut_type_plot.rds'), path('*_driver_oncoprint.rds'), emit: data_stats_rds
    tuple val(meta), path('*_mutations_report.pdf'), emit: mutation_report_pdf
    tuple val(meta), path('*_driver_oncoprint.pdf'), emit: driver_oncoprint_plot_pdf
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    template "main_script.R"

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    touch ${prefix}_vaf_chr_plot.rds
    touch ${prefix}_chr_mut.rds
    touch ${prefix}_consequence_mut_plot.rds
    touch ${prefix}_mut_type_plot.rds
    touch ${prefix}_driver_oncoprint.rds
    touch ${prefix}_mutations_report.pdf
    touch ${prefix}_driver_oncoprint.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dplyr: \$(Rscript -e "library(dplyr); cat(as.character(packageVersion('dplyr')))")
        tidyr: \$(Rscript -e "library(tidyr); cat(as.character(packageVersion('tidyr')))")
        ggplot2: \$(Rscript -e "library(ggplot2); cat(as.character(packageVersion('ggplot2')))")
        complexheatmap: \$(Rscript -e "library(ComplexHeatmap); cat(as.character(packageVersion('ComplexHeatmap')))")
        patchwork: \$(Rscript -e "library(patchwork); cat(as.character(packageVersion('patchwork')))")
    END_VERSIONS
    """
}
