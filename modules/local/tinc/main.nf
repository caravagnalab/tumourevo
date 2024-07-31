process TINC {

    tag "$meta.id"
    container 'docker://vvvirgy/tinc:v2' // define the running container
    // conda "${moduleDir}/environment.yml"
  
    input:
   
    tuple val(meta), path(cna_RDS)
    tuple val(meta), path(snv_RDS)
    // tuple val(datasetID), val(patientID), val(sampleID), path(snv_RDS)
  
  output:
    tuple val(meta), path("*_fit.rds"), emit: rds 
    tuple val(meta), path("*_plot.rds"), emit: plot_rds
    tuple val(meta), path("*.pdf"), emit: plot_pdf
    

  script:

    def args                                = task.ext.args 
    def prefix                              = task.ext.prefix                                       ?: "${meta.id}"
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

    all_mutations = readRDS("$snv_RDS")
    samples = names(all_mutations)

    tumor_sample = "$meta.tumour_sample"
    normal_sample = "$meta.normal_sample"
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
    
    TINC_fit = TINC::autofit(input = input_mut, 
                    cna = CNAs, 
                    VAF_range_tumour = eval(parse(text="$vaf_range_tumour")),
                    cutoff_miscalled_clonal = eval(parse(text="$cutoff_miscalled_clonal")),
                    cutoff_lv_assignment = eval(parse(text="$cutoff_lv_assignment")),
                    N = eval(parse(text="$N")),
                    FAST = eval(parse(text="$fast"))
                    )
                    
    tinc_plot = plot(TINC_fit)
    
    saveRDS(filename = paste0("${prefix}", "_plot.rds"), object = tinc_plot)
    ggplot2::ggsave(plot = tinc_plot, filename = paste("${prefix}", "_plot.pdf"), width = 210, height = 297, units="mm", dpi = 200)
    saveRDS(filename = paste0("${prefix}", "_fit.rds"), object = TINC_fit)

    """
}
