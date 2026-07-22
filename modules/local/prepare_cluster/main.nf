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

    convert_ccf_to_vaf <- function(ccf, purity, major_cn, minor_cn, multiplicity, cn_normal = 2) {
      cn_tumor <- major_cn + minor_cn
      vaf <- (ccf * multiplicity * purity) / ((1 - purity) * cn_normal + purity * cn_tumor)
      pmin(pmax(vaf, eps), 1 - eps)
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

	clusters = clonal_clusters %>% pull(cluster) %>% unique()
        data_table <- readr::read_tsv("$data") %>%
            dplyr::select(chr, from, to, ref, alt, is_driver, driver_label, QC_PASS, blacklisted, karyotype, DP, NV, VAF, Indiv, CCF, mutation_multiplicity, purity) %>%
            dplyr::mutate(chr = sub("chr", "", chr)) %>%
            dplyr::mutate(mutation_id = paste("${meta.patient}", chr, from, alt, sep = ':'))

        driver <- data_table %>% dplyr::filter(is_driver == T & karyotype == '1:1')
        driver_cn <- data_table %>% dplyr::filter(is_driver == T & karyotype != '1:1')


        mut_data = tool_table\$data %>%
          bind_cols(tool_table\$labels) %>%
          tidyr::pivot_longer(
            cols = c(matches("^NV"), matches("^DP"), matches("^VAF")),
            names_to = c(".value", "sample"),
            names_sep = "\\\\."
          ) %>%
          select(chr, from, ref, alt, gene, driver, NV, DP, VAF, sample, cluster.Binomial) %>%
          mutate(mutation_id = paste("${meta.patient}", sub("^chr", "", chr), from, alt, sep = ':'))

        to_assign_muts <- driver %>% filter(!mutation_id %in% mut_data\$mutation_id)
        to_assign_muts_cn <- driver_cn %>% filter(!mutation_id %in% mut_data\$mutation_id)

        if (nrow(to_assign_muts)>0){
          #theta = tool_table\$theta_k %>% as.data.frame() %>% tibble::add_column(sample = rownames(.))
          theta = tool_table\$theta_k[,clusters] %>% as.data.frame() %>% tibble::add_column(sample = rownames(.))
	  samples = theta\$sample %>% unique()
          clusters = colnames(theta)[1:(ncol(theta)-1)]
          eps = 1e-9

          to_assign_muts = to_assign_muts %>% dplyr::rename(sample = Indiv)
          mut_data_sample = to_assign_muts %>% filter(sample %in% samples) #filter(sample == !!samples)

          mut_cluster_assignments <- mut_data_sample %>%
            rowwise() %>%
            group_by(chr, from, ref, alt, is_driver, driver_label) %>%
            summarize(
              cluster = {
                log_liks <- setNames(
                  vapply(clusters, function(k) {
                    mu <- pmax(pmin(theta[as.character(samples), k], 1 - eps), eps)
                    sum(dbinom(x = NV, size = DP, prob = mu, log = TRUE))
                  }, numeric(1)),
                  clusters
                )
                names(which.max(log_liks))
              }
            ) %>%
            ungroup() %>%
            left_join(clonal_clusters %>% select(cluster, is_clonal) %>% mutate(cluster = as.character(cluster)) %>% distinct()) %>%
            select(chr, from, ref, alt, is_driver, driver_label, cluster, is_clonal) %>%
            mutate(chr = paste0('chr', chr))
        } else {
          mut_cluster_assignments <- tibble()
        }

        if (nrow(to_assign_muts_cn)>0){
          clonal = unique(clonal_clusters %>% filter(is_clonal == T) %>% pull(cluster))
          theta = tool_table\$theta_k %>% as.data.frame() %>% tibble::add_column(sample = rownames(.))
          samples = theta\$sample %>% unique()
          theta_clonal = theta[,clonal]
          names(theta_clonal) = rownames(theta)
          purity = to_assign_muts_cn\$purity %>% unique()
          ccf_clonal = (theta_clonal * ((1+1-2)*purity+2)) / (1 * purity)
          eps = 1e-9

          to_assign_muts_cn = to_assign_muts_cn %>% dplyr::rename(sample = Indiv)
          mut_data_sample_cn = to_assign_muts_cn %>%
            filter(sample %in% samples) %>%
            tidyr::separate(col = karyotype, into = c('major', 'minor'), sep = ':', remove = F, convert = T)

          mut_cluster_assignments_cn <- mut_data_sample_cn %>%
            rowwise() %>%
            group_by(chr, from, ref, alt, is_driver, driver_label) %>%
            summarize(
              log_lik_clonal = {
                mu_tmp <- (mutation_multiplicity * purity * ccf_clonal[sample])/(2*(1-purity) + purity*(major+minor))
                mu <- pmax(pmin(mu_tmp, 1 - eps), eps)
                sum(dbinom(x = NV, size = DP, prob = mu, log = TRUE))
              }
            ) %>%
            ungroup() %>%
            mutate(cluster = ifelse(log_lik_clonal > -3, clonal, NA)) %>%
            mutate(is_clonal = ifelse(log_lik_clonal > -3, T, F)) %>%
            select(chr, from, ref, alt, is_driver, driver_label, cluster, is_clonal) %>%
            mutate(chr = paste0('chr', chr)) %>%
            filter(is_clonal) %>%
            mutate(cluster = as.character(cluster))
        } else {
          mut_cluster_assignments_cn <- tibble()
        }


        sigprofiler_table <- tool_table\$data %>%
            dplyr::mutate(cluster=tool_table\$labels\$cluster.Binomial) %>%
            dplyr::left_join(clonal_clusters) %>%
            dplyr::rename(driver_label=gene, is_driver=driver) %>%
            dplyr::select(chr, from, ref, alt, cluster, driver_label, is_driver, is_clonal) %>%
            dplyr::bind_rows(mut_cluster_assignments) %>%
            dplyr::bind_rows(mut_cluster_assignments_cn) %>%
            dplyr::distinct() %>%
            tidyr::separate(col = chr, sep = 'chr', into = c('tmp', 'chr')) %>%
            dplyr::select(-tmp) %>%
            dplyr::mutate(Project="${meta.id}", Genome="${params.genome}", mut_type='SNP', Type='SOMATIC', ID=cluster, Sample=cluster) %>%
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
            dplyr::mutate(Project = "${meta.id}", Genome = "${params.genome}", mut_type = 'SNP', Type = 'SOMATIC', ID = cluster, Sample = cluster) %>%
            dplyr::rename(chrom = chr, pos_start = from) %>%
            dplyr::rowwise() %>%
            dplyr::mutate(pos_end = pos_start + abs(stringr::str_count(ref) - stringr::str_count(alt))) %>%
            dplyr::select(Project, Sample, ID, Genome, mut_type, chrom, pos_start, pos_end, ref,alt, Type,  driver_label, is_driver, is_clonal) %>%
            dplyr::filter(ref != alt) %>%
            dplyr::distinct()

    } else if (grepl(pattern = 'best_fit.txt', x = "$fit")){
        tool <- 'pyclonevi'

        best_fit_table = readr::read_tsv("$fit") %>%
            dplyr::rename(ccf=cellular_prevalence, cluster=cluster_id) %>%
            group_by(cluster) %>%
            filter(!all(ccf == 0)) %>%
            ungroup()

        clonal_clusters = best_fit_table %>%
            dplyr::mutate(is_clonal=ifelse(cluster==get_clonal_cluster(best_fit_table),
                                           TRUE, FALSE)) %>%
            dplyr::select(mutation_id, cluster, is_clonal) %>% unique()

        cluster_table <- best_fit_table %>%
            dplyr::select(mutation_id, cluster) %>%
            dplyr::distinct()

        cluster_ccf <- best_fit_table %>% select(sample_id, mutation_id, cluster, ccf) %>% distinct()


        data_table <- readr::read_tsv("$data") %>%
            dplyr::select(chr, from, to, ref, alt, is_driver, driver_label, QC_PASS, blacklisted, karyotype, DP, NV, VAF, Indiv, purity) %>%
            dplyr::mutate(chr = sub("chr", "", chr)) %>%
            dplyr::mutate(mutation_id = paste("$meta.patient", chr, from, alt, sep = ':'))

        driver <- data_table %>% dplyr::filter(is_driver == T)

        to_assign_muts <- driver %>% filter(!mutation_id %in% cluster_table\$mutation_id)


        if (nrow(to_assign_muts)>0){

          samples = cluster_ccf\$sample_id %>% unique()
          clusters = cluster_ccf\$cluster %>% unique()
          eps = 1e-9

          to_assign_muts = to_assign_muts %>% dplyr::rename(sample = Indiv)
          mut_data_sample = to_assign_muts %>% filter(sample == !!samples)

          pi = data_table\$purity %>% unique()

          mut_cluster_assignments <- mut_data_sample %>%
            tidyr::separate(karyotype, into = c('Major', 'minor'), sep = ':', remove = F, convert = T) %>%
            rowwise() %>%
            group_by(chr, from, ref, alt, is_driver, driver_label) %>%
            summarize(
              cluster = {
                log_liks <- setNames(
                  vapply(clusters, function(k) {
                    ccf_cluster <- cluster_ccf %>% filter(cluster == k) %>% pull(ccf) %>% unique()
                    mu <- convert_ccf_to_vaf(ccf_cluster, purity = pi, major_cn = Major, minor_cn = minor, multiplicity = 1)
                    sum(dbinom(x = NV, size = DP, prob = mu, log = TRUE))
                  }, numeric(1)),
                  clusters
                )
                names(which.max(log_liks))
              }
            ) %>%
            ungroup() %>%
            left_join(clonal_clusters %>% select(cluster, is_clonal) %>% mutate(cluster = as.character(cluster)) %>% distinct()) %>%
            select(chr, from, ref, alt, is_driver, driver_label, cluster, is_clonal) %>%
            mutate(cluster = as.numeric(cluster))
        } else {
          mut_cluster_assignments <- tibble()
        }


        sigprofiler_table <- cluster_table %>%
            dplyr::left_join(data_table) %>%
            dplyr::left_join(clonal_clusters) %>%
            dplyr::select(chr, from, ref, alt, is_driver, driver_label, cluster, is_clonal) %>%
            dplyr::bind_rows(mut_cluster_assignments) %>%
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
