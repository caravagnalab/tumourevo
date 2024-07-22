process TINC {

    tag "$meta.id"
    // label 'process_medium' // magari ha senso iniziare ad aggiungerlo?
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/gawk:5.1.0' :
    //     'biocontainers/gawk:5.1.0' }"
    container 'docker://vvvirgy/tinc:v2' // define the running container
    // conda "${moduleDir}/environment.yml"
    // publishDir params.publish_dir, mode: 'copy'

    input:
   
    tuple val(meta), path(cna_RDS), path(snv_RDS)
    // tuple val(datasetID), val(patientID), val(sampleID), path(snv_RDS)
  
  output:
    tuple val(meta), path("*_fit.rds"), emit: rds 
    tuple val(meta), path("*_plot.rds"), emit: plot_rds
    tuple val(meta), path("*.pdf"), emit: plot_pdf
    // tuple val(datasetID), val(patientID), val(sampleID), path("QC/TINC/$datasetID/$patientID/$sampleID/TINC/*tinc_plot.rds"), emit: plot_rds
    // tuple val(datasetID), val(patientID), val(sampleID), path("QC/TINC/$datasetID/$patientID/$sampleID/TINC/*tinc_fit.rds"), emit: fit_rds
    // tuple val(datasetID), val(patientID), val(sampleID), path("QC/TINC/$datasetID/$patientID/$sampleID/TINC/*plot.pdf"), emit: plot_pdf

  script:

    def args                                = task.ext.args 
    def prefix                              = task.ext.prefix                                       ?: "${meta.id}"
    def matching_strategy                   = args!='' && args.matching_strategy                    ?  "$args.matching_strategy" : ""
    def vaf_range_tumour                    = args!='' && args.vaf_range_tumour                     ?  "$args.vaf_range_tumour"  : ""
    def cutoff_miscalled_clonal             = args!='' && args.cutoff_miscalled_clonal              ?  "$args.cutoff_miscalled_clonal" : ""
    def cutoff_lv_assignment                = args!='' && args.cutoff_lv_assignment                 ?  "$args.cutoff_lv_assignment"  : ""
    def N                                   = args!='' && args.N                                    ?  "$args.N"  : ""
    def fast                                = args!='' && args.fast                                 ?  "$args.fast" : ""
    

    """
    #!/usr/bin/env Rscript

    require(CNAqc)
    require(tidyverse)
    require(TINC)

    # res_dir = paste0("QC/TINC/", "$datasetID", "/", "$patientID", "/", "$sampleID") --> capire se la directory viene creata prima o meno
    # dir.create(res_dir, recursive = TRUE)

    all_mutations = readRDS("$snv_RDS")
    samples = names(all_mutations)

    tumor_sample = "$meta.sampleID"
    normal_sample = setdiff(samples, tumor_sample)
    tumor_mutations = all_mutations[[tumor_sample]]\$mutations %>% 
      select(chr, from, to, ref, alt, NV, DP, NR, VAF) %>% 
      rename(t_alt_count = NV, t_ref_count = NR, t_tot_count = DP, t_vaf = VAF)

    normal_mutations = all_mutations[[normal_sample]]\$mutations %>% 
      select(chr, from, to, ref, alt, NV, DP, NR,VAF) %>% 
      rename(n_alt_count = NV, n_ref_count = NR, n_tot_count = DP, n_vaf = VAF)

    input_mut <- dplyr::full_join(tumor_mut, normal_mut, by = c("chr", "from", "ref", "alt")) %>% 
        mutate(VAF.x = case_when(is.na(VAF.x) ~ 0, .default = VAF.x)) %>%
        mutate(VAF.y = case_when(is.na(VAF.y) ~ 0, .default = VAF.y)) %>%
        mutate(VAF.x = as.numeric(VAF.x), VAF.y = as.numeric(VAF.y))

    CNAs = readRDS("$cna_RDS")\$segments
    
    TINC_fit = TINC::autofit(input_mut, cna = CNAs, FAST = FALSE)
    tinc_plot = plot(TINC_fit)
    
    saveRDS(filename = paste0("${prefix}", "_plot.rds"), object = tinc_plot)
    ggplot2::ggsave(plot = tinc_plot, filename = paste("${prefix}", "_plot.pdf"), width = 210, height = 297, units="mm", dpi = 200)
    saveRDS(filename = paste0("${prefix}", "_fit.rds"), object = TINC_fit)

    """
}
