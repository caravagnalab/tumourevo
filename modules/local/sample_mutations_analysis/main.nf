process SAMPLE_MUTATIONS_ANALYSIS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f9/f999ab0c279823df3a494cf13dda7fdbd4665f1f62563caafa7e5f67811859b7/data':
        'community.wave.seqera.io/library/bioconductor-complexheatmap_bioconductor-maftools_r-dplyr_r-ggplot2_pruned:0431e79bdf418942' }"
        
    input:
    tuple val(meta), path(snv_rds), val(tumour_sample)
    tuple val(meta2), path(maf)

    output:
    tuple val(meta), path("*_tmb.rds"), emit: tmb_rds 
    tuple val(meta), path("*_vaf_chr_plot.rds"), path('*_chr_mut.rds'), path('*_driver_oncoprint.rds'), emit: data_stats_rds
    tuple val(meta), path('*_mutations_report.pdf'), emit: mutation_report_pdf
    tuple val(meta), path('*_mutations_per_chr.pdf'), emit: mutations_per_chr
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
    touch ${prefix}_driver_oncoprint.rds
    touch ${prefix}_mutations_report.pdf
    touch ${prefix}_mutations_per_chr.pdf
    touch ${prefix}_driver_oncoprint.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dplyr: \$(Rscript -e "library(dplyr); cat(as.character(packageVersion('dplyr')))")
        tidyr: \$(Rscript -e "library(tidyr); cat(as.character(packageVersion('tidyr')))")
        ggplot2: \$(Rscript -e "library(ggplot2); cat(as.character(packageVersion('ggplot2')))")
        complexheatmap: \$(Rscript -e "library(ComplexHeatmap); cat(as.character(packageVersion('ComplexHeatmap')))")
        patchwork: \$(Rscript -e "library(patchwork); cat(as.character(packageVersion('patchwork')))")
        maftools: \$(Rscript -e "library(maftools); cat(as.character(packageVersion('maftools')))")
    END_VERSIONS
    """
}
