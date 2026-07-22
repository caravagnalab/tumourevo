process COHORT_MUTATIONS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e0/e057edad2ef577776755a1397e4831ebe1a375617ad08dab1a0fb2044a3c06c6/data':
        'community.wave.seqera.io/library/bioconductor-complexheatmap_bioconductor-maftools_r-cnaqc_r-dplyr_pruned:ab37afa63a496cf3' }"
        
    input:
    tuple val(meta), path(join_cnaqc), val(cnaqc_patients), path(tmb_rds), val(tmb_patients), path(maf)

    output:
    tuple val(meta), path("*_oncoprint_tmb.pdf"), path("*_tmb.pdf"), emit: mutation_report
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
    touch ${prefix}_oncoprint.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dplyr: \$(Rscript -e "library(dplyr); cat(as.character(packageVersion('dplyr')))")
        tidyr: \$(Rscript -e "library(tidyr); cat(as.character(packageVersion('tidyr')))")
        ggplot2: \$(Rscript -e "library(ggplot2); cat(as.character(packageVersion('ggplot2')))")
        complexheatmap: \$(Rscript -e "library(ComplexHeatmap); cat(as.character(packageVersion('ComplexHeatmap')))")
        cnaqc: \$(Rscript -e "library(CNAqc); cat(as.character(packageVersion('CNAqc')))")
        maftools: \$(Rscript -e "library(maftools); cat(as.character(packageVersion('maftools')))")
        cowplot: \$(Rscript -e "library(cowplot); cat(as.character(packageVersion('cowplot')))")
    END_VERSIONS
    """
}
