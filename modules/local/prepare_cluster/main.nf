process PREPARE_CLUSTER {
    tag "$meta.id"
    label "process_single"
    label "error_retry"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/91/913112a2d9295d35fe97874caf5f500df5a98d6fef1cb5861fd64caa0223a047/data':
        'community.wave.seqera.io/library/r-cnaqc_r-cli_r-dplyr_r-readr_pruned:0fc82bfd06afe6dc' }"

    input:
    tuple val(meta), path(fit), path(data)

    output:
    //tuple val(meta), path("*_all_positions.rds"), emit: all_pos
    path "versions.yml",                          emit: versions

    script:
    def args   = task.ext.args   ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    #!/usr/bin/env Rscript
    library(dplyr)

    if (grepl(pattern = 'viber', x = "$fit")){
        print('is viber')
        tool_table <- readRDS("$fit")

        sigprofiler_table <- tool_table\$data %>%
            dplyr::select(chr, from, ref, alt) %>%
            dplyr::distinct() %>%
            tidyr::separate(col = chr, sep = 'chr', into = c('tmp', 'chr')) %>%
            dplyr::select(-tmp) %>%
            dplyr::bind_cols(tool_table\$labels) %>%
            dplyr::rename(cluster = cluster.Binomial) %>%
            dplyr::mutate(Project = "$meta.id", Genome = 'GRCh38', mut_type = 'SNP', Type = 'SOMATIC', ID = cluster, Sample = cluster) %>%
            dplyr::rename(chrom = chr, pos_start = from) %>%
            dplyr::rowwise() %>%
            dplyr::mutate(pos_end = pos_start + abs(stringr::str_count(ref) - stringr::str_count(alt))) %>%
            dplyr::select(Project, Sample, ID, Genome, mut_type, chrom, pos_start, pos_end, ref,alt, Type) %>%
            dplyr::filter(ref != alt) %>%
            dplyr::distinct()


    } else if (grepl(pattern = 'mobster', x = "$fit")){
        tool_table <- readRDS("$fit")

        sigprofiler_table <- tool_table\$data %>%
            dplyr::select(chr, from, ref, alt, cluster) %>%
            dplyr::distinct() %>%
            tidyr::separate(col = chr, sep = 'chr', into = c('tmp', 'chr')) %>%
            dplyr::select(-tmp)  %>%
            dplyr::mutate(Project = "$meta.id", Genome = 'GRCh38', mut_type = 'SNP', Type = 'SOMATIC', ID = cluster, Sample = cluster) %>%
            dplyr::rename(chrom = chr, pos_start = from) %>%
            dplyr::rowwise() %>%
            dplyr::mutate(pos_end = pos_start + abs(stringr::str_count(ref) - stringr::str_count(alt))) %>%
            dplyr::select(Project, Sample, ID, Genome, mut_type, chrom, pos_start, pos_end, ref,alt, Type) %>%
            dplyr::filter(ref != alt) %>%
            dplyr::distinct()

    } else if (grepl(pattern = 'best_fit.txt', x = "$fit")){
        cluster_table <- readr::read_tsv("$fit") %>%
            dplyr::select(mutation_id, cluster_id) %>%
            dplyr::distinct()
        print(cluster_table)

        data_table <- readr::read_tsv("$data") %>%
            dplyr::select(chr, from, to, ref, alt) %>%
            dplyr::mutate(chr = sub("chr", "", chr)) %>%
            dplyr::mutate(mutation_id = paste("$meta.patient", chr, from, alt, sep = ':'))
        print(data_table)

        sigprofiler_table <- cluster_table %>%
            dplyr::left_join(data_table) %>%
            dplyr::select(-mutation_id) %>%
            dplyr::rename(cluster = cluster_id) %>%
            dplyr::mutate(Project = "$meta.id", Genome = 'GRCh38', mut_type = 'SNP', Type = 'SOMATIC', ID = cluster, Sample = cluster) %>%
            dplyr::rename(chrom = chr, pos_start = from) %>%
            dplyr::rowwise() %>%
            dplyr::mutate(pos_end = pos_start + abs(stringr::str_count(ref) - stringr::str_count(alt))) %>%
            dplyr::select(Project, Sample, ID, Genome, mut_type, chrom, pos_start, pos_end, ref,alt, Type) %>%
            dplyr::filter(ref != alt) %>%
            dplyr::distinct()

    }


    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_all_positions.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dplyr: \$(Rscript -e "cat(as.character(packageVersion('dplyr')))")
    END_VERSIONS
    """
}
