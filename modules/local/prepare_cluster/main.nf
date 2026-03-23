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
    tuple val(meta), path("*.txt"), emit: signature_table
    path "versions.yml",            emit: versions

    script:
    def args   = task.ext.args   ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    #!/usr/bin/env Rscript
    library(dplyr)

    get_clonal_cluster = function(df) {
        theta_long = df %>%
            dplyr::group_by(sample_id, cluster) %>%
            dplyr::summarize(ccf=mean(ccf, na.rm=TRUE), .groups="drop")
        theta = theta_long %>%
            tidyr::pivot_wider(names_from=cluster, values_from=ccf) %>%
            dplyr::select(-sample_id)
        max_colnames = apply(theta, 1, function(row) {
            names(row)[which(row == max(row))]
        }) # Extract all clusters which have max ccf for each sample (because in one sample there can be more than one cluster with ccf == 1)
        names(which.max(table(unlist(max_colnames)))) # extract the cluster which appear more frequently (i.e. possibly in all the samples)
    }

    if (grepl(pattern = 'viber', x = "$fit")){
        tool <- 'viber'
        tool_table <- readRDS("$fit")

        clonal_clusters = tool_table\$data %>% 
            dplyr::mutate(cluster=tool_table\$labels\$cluster.Binomial) %>%
            tidyr::pivot_longer(cols=starts_with("VAF"), names_to="sample_id",
                                values_to="ccf", names_prefix="VAF.")
        clonal_clusters = clonal_clusters %>%
            dplyr::mutate(is_clonal=ifelse(cluster==get_clonal_cluster(clonal_clusters),
                                           TRUE, FALSE)) %>%
            dplyr::select(chr, from, ref, alt, cluster, is_clonal) %>% unique()

        sigprofiler_table <- tool_table\$data %>%
            dplyr::mutate(cluster=tool_table\$labels\$cluster.Binomial) %>%
            dplyr::left_join(clonal_clusters) %>%
            dplyr::rename(driver_label=gene, is_driver=driver) %>%
            dplyr::select(chr, from, ref, alt, cluster, driver_label, is_driver, is_clonal) %>%
            dplyr::distinct() %>%
            tidyr::separate(col = chr, sep = 'chr', into = c('tmp', 'chr')) %>%
            dplyr::select(-tmp) %>%
            dplyr::mutate(Project="$meta.id", Genome="${params.genome}", mut_type='SNP', Type='SOMATIC', ID=cluster, Sample=cluster) %>%
            dplyr::rename(chrom = chr, pos_start = from) %>%
            dplyr::rowwise() %>%
            dplyr::mutate(pos_end = pos_start + abs(stringr::str_count(ref) - stringr::str_count(alt))) %>%
            dplyr::select(Project, Sample, ID, Genome, mut_type, chrom, pos_start, pos_end, ref, alt, Type, driver_label, is_driver, is_clonal) %>%
            dplyr::filter(ref != alt) %>%
            dplyr::distinct()

    } else if (grepl(pattern = 'mobster', x = "$fit")){
        tool <- 'mobster'
        tool_table <- readRDS("$fit")

        clonal_clusters = tool_table\$data %>%
            dplyr::mutate(is_clonal=ifelse(cluster==get_clonal_cluster(tool_table\$data %>% dplyr::rename(ccf=VAF)),
                                           TRUE, FALSE)) %>%
            dplyr::select(chr, from, ref, alt, cluster, is_clonal) %>% unique()

        sigprofiler_table <- tool_table\$data %>%
            dplyr::select(chr, from, ref, alt, cluster, driver_label, is_driver) %>%
            dplyr::left_join(clonal_clusters) %>%
            dplyr::distinct() %>%
            tidyr::separate(col = chr, sep = 'chr', into = c('tmp', 'chr')) %>%
            dplyr::select(-tmp)  %>%
            dplyr::mutate(Project = "$meta.id", Genome = "${params.genome}", mut_type = 'SNP', Type = 'SOMATIC', ID = cluster, Sample = cluster) %>%
            dplyr::rename(chrom = chr, pos_start = from) %>%
            dplyr::rowwise() %>%
            dplyr::mutate(pos_end = pos_start + abs(stringr::str_count(ref) - stringr::str_count(alt))) %>%
            dplyr::select(Project, Sample, ID, Genome, mut_type, chrom, pos_start, pos_end, ref,alt, Type,  driver_label, is_driver, is_clonal) %>%
            dplyr::filter(ref != alt) %>%
            dplyr::distinct()

    } else if (grepl(pattern = 'best_fit.txt', x = "$fit")){
        tool <- 'pyclonevi'

        best_fit_table = readr::read_tsv("$fit") %>%
            dplyr::rename(ccf=cellular_prevalence, cluster=cluster_id)

        clonal_clusters = best_fit_table %>%
            dplyr::mutate(is_clonal=ifelse(cluster==get_clonal_cluster(best_fit_table),
                                           TRUE, FALSE)) %>%
            dplyr::select(mutation_id, cluster, is_clonal) %>% unique()

        cluster_table <- best_fit_table %>%
            dplyr::select(mutation_id, cluster) %>%
            dplyr::distinct()

        data_table <- readr::read_tsv("$data") %>%
            dplyr::select(chr, from, to, ref, alt, is_driver, driver_label) %>%
            dplyr::mutate(chr = sub("chr", "", chr)) %>%
            dplyr::mutate(mutation_id = paste("$meta.patient", chr, from, alt, sep = ':'))

        sigprofiler_table <- cluster_table %>%
            dplyr::left_join(data_table) %>%
            dplyr::left_join(clonal_clusters) %>%
            dplyr::select(-mutation_id) %>%
            dplyr::mutate(cluster = paste0('C', cluster)) %>%
            dplyr::mutate(Project = "$meta.id", Genome = "${params.genome}", mut_type = 'SNP', Type = 'SOMATIC', ID = cluster, Sample = cluster) %>%
            dplyr::rename(chrom = chr, pos_start = from) %>%
            dplyr::rowwise() %>%
            dplyr::mutate(pos_end = pos_start + abs(stringr::str_count(ref) - stringr::str_count(alt))) %>%
            dplyr::select(Project, Sample, ID, Genome, mut_type, chrom, pos_start, pos_end, ref,alt, Type,  driver_label, is_driver, is_clonal) %>%
            dplyr::filter(ref != alt) %>%
            dplyr::distinct()
    }

    write.table(sigprofiler_table, file = paste0("$prefix", "_", tool, ".txt"), quote = F, sep = '\\t', row.names = F)

    # version export
    f <- file("versions.yml","w")
    dplyr_version <- sessionInfo()\$otherPkgs\$dplyr\$Version
    tidyr_version <- sessionInfo()\$otherPkgs\$tidyr\$Version
    vcfR_version <- sessionInfo()\$otherPkgs\$vcfR\$Version
    writeLines(paste0('"', "$task.process", '"', ":"), f)
    writeLines(paste("    dplyr:", dplyr_version), f)
    writeLines(paste("    tidyr:", tidyr_version), f)
    close(f)
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dplyr: \$(Rscript -e "cat(as.character(packageVersion('dplyr')))")
        tidyr: \$(Rscript -e "cat(as.character(packageVersion('tidyr')))")
    END_VERSIONS
    """
}
