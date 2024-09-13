process TINC {

    tag "$meta.id"
    container 'docker://vvvirgy/tinc:v2' // define the running container
    // conda "${moduleDir}/environment.yml"
  
    input:
   
    tuple val(meta), path(cna_RDS), path(snv_RDS)
  
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
      filter(!is.na(DP)) %>%
      rename(t_alt_count = NV, t_ref_count = NR, t_tot_count = DP, t_vaf = VAF)

    normal_mutations = all_mutations[[normal_sample]]\$mutations %>% 
      select(chr, from, to, ref, alt, NV, DP, NR, VAF) %>% 
      filter(!is.na(DP)) %>%
      rename(n_alt_count = NV, n_ref_count = NR, n_tot_count = DP, n_vaf = VAF)

    input_mut = dplyr::full_join(tumor_mutations, normal_mutations, by = c("chr", "from", "to", "ref", "alt")) %>%
        mutate(t_vaf = case_when(is.na(t_vaf) ~ 1e-5, .default = t_vaf)) %>%
        mutate(n_vaf = case_when(is.na(n_vaf) ~ 1e-5, .default = n_vaf)) %>%
        mutate(t_vaf = as.numeric(t_vaf), n_vaf = as.numeric(n_vaf)) %>%
        mutate(t_alt_count = case_when(is.na(t_alt_count) ~ 0, .default = t_alt_count)) %>%
        mutate(t_ref_count = case_when(is.na(t_ref_count) ~ 0, .default = t_ref_count)) %>%
        mutate(n_alt_count = case_when(is.na(n_alt_count) ~ 0, .default = n_alt_count)) %>%
        mutate(n_ref_count = case_when(is.na(n_ref_count) ~ 0, .default = n_ref_count)) %>%
        mutate(t_alt_count = as.numeric(t_alt_count), t_ref_count = as.numeric(t_ref_count), n_tot_count = as.numeric(n_tot_count), n_ref_count = as.numeric(n_ref_count))

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
    
    saveRDS(file = paste0("${prefix}", "_plot.rds"), object = tinc_plot)
    ggplot2::ggsave(plot = tinc_plot, filename = paste("${prefix}", "_plot.pdf"), width = 210, height = 297, units="mm", dpi = 200)
    saveRDS(file = paste0("${prefix}", "_fit.rds"), object = TINC_fit)

    """
}
