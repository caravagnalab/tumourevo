process CNAQC {
  tag "$meta.id"
  container = 'docker://lvaleriani/cnaqc:dev1'

  input:
    
    tuple val(meta), path(cna_rds)
    tuple val(meta), path(snv_rds)
  
  output:

    tuple val(meta), path("*.rds"), emit: qc_rds
    tuple val(meta), path("*.rds"), path("*.rds"), emit: plot_rds
    tuple val(meta), path("*.pdf"), emit: plot_pdf_data
    tuple val(meta), path("*.pdf"), emit: plot_pdf_qc

  script:

    def args                                = task.ext.args    
    def prefix                              = task.ext.prefix                                       ?: "${meta.id}"  
    def matching_strategy                   = args!='' && args.matching_strategy                    ?  "$args.matching_strategy" : ""
    def karyotypes                          = args!='' && args.karyotypes                           ?  "$args.karyotypes" : ""
    def min_karyotype_size                  = args!='' && args.min_karyotype_size                   ?  "$args.min_karyotype_size" : ""
    def min_absolute_karyotype_mutations    = args!='' && args.min_absolute_karyotype_mutations     ?  "$args.min_absolute_karyotype_mutations" : ""
    def p_binsize_peaks                     = args!='' && args.p_binsize_peaks                      ?  "$args.p_binsize_peaks" : ""
    def matching_epsilon                    = args!='' && args.matching_epsilon                     ?  "$args.matching_epsilon" : ""
    def purity_error                        = args!='' && args.purity_error                         ?  "$args.purity_error" : ""
    def vaf_tolerance                       = args!='' && args.vaf_tolerance                        ?  "$args.vaf_tolerance" : ""
    def n_bootstrap                         = args!='' && args.n_bootstrap                          ?  "$args.n_bootstrap" : ""
    def kernel_adjust                       = args!='' && args.kernel_adjust                        ?  "$args.kernel_adjust" : ""
    def kde                                 = args!='' && args.kde                                  ?  "$args.KDE" : ""
    def starting_state_subclonal_evolution  = args!='' && args.starting_state_subclonal_evolution   ?  "$args.starting_state_subclonal_evolution" : ""
    def cluster_subclonal_CCF               = args!='' && args.cluster_subclonal_CCF                ?  "$args.cluster_subclonal_CCF" : ""

    def muts_per_karyotype                  = args!='' && args.muts_per_karyotype                   ?  "$args.muts_per_karyotype" : ""
    def cutoff_QC_PASS                      = args!='' && args.cutoff_QC_PASS                       ?  "$args.cutoff_QC_PASS" : ""
    def method                              = args!='' && args.method                               ?  "$args.method" : ""

    def plot_cn                             = args!='' && args.plot_cn                               ?  "$args.plot_cn" : ""
    """
    #!/usr/bin/env Rscript

    library(tidyverse)
    library(CNAqc)
    
    SNV = readRDS("$snv_rds")
    SNV = SNV[["$meta.tumour_sample"]]
    SNV = SNV\$mutations 
    SNV = SNV %>% dplyr::mutate(mutation_id = paste(chr,from,to,ref,alt,sep = ':'))

    CNA = readRDS("$cna_rds")
    x = CNAqc::init(
        mutations = SNV,
        cna = CNA\$segments,
        purity = CNA\$purity ,
        ref = "$params.genome")

    x = CNAqc::analyze_peaks(x, 
      matching_strategy = "$matching_strategy",
      karyotypes =  eval(parse(text="$karyotypes")),
      min_karyotype_size = as.numeric("$min_karyotype_size"),
      min_absolute_karyotype_mutations = as.numeric("$min_absolute_karyotype_mutations"),
      p_binsize_peaks = as.numeric("$p_binsize_peaks"),
      matching_epsilon = eval(parse(text="$matching_epsilon")),
      purity_error = as.numeric("$purity_error"),
      VAF_tolerance = as.numeric("$vaf_tolerance"),
      n_bootstrap = as.numeric("$n_bootstrap"),
      kernel_adjust = as.numeric("$kernel_adjust"),
      KDE = eval(parse(text = "$kde")),
      starting_state_subclonal_evolution = "$starting_state_subclonal_evolution",
      cluster_subclonal_CCF = as.logical("$cluster_subclonal_CCF"),
      min_VAF = 0
      )

    x = CNAqc::compute_CCF(
      x,
      karyotypes = eval(parse(text="$karyotypes")),
      muts_per_karyotype = as.numeric("$muts_per_karyotype"),
      cutoff_QC_PASS = as.numeric("$cutoff_QC_PASS"),
      method = "$method"
    )

    pl = ggpubr::ggarrange(
      CNAqc::plot_data_histogram(x, which = 'VAF', karyotypes = eval(parse(text="$karyotypes"))),
      CNAqc::plot_data_histogram(x, which = 'DP', karyotypes = eval(parse(text="$karyotypes"))),
      CNAqc::plot_data_histogram(x, which = 'NV', karyotypes = eval(parse(text="$karyotypes"))),
      CNAqc::plot_data_histogram(x, which = 'CCF', karyotypes = eval(parse(text="$karyotypes"))),
      ncol = 2,
      nrow = 2
    )


    pl_exp = patchwork::wrap_plots(
      CNAqc::plot_gw_counts(x),
      CNAqc::plot_gw_vaf(x, N = 10000),
      CNAqc::plot_gw_depth(x, N = 10000),
      CNAqc::plot_segments(x),
      pl, 
      nrow = 5, 
      heights = c(.5,.5,.5,1,5)) + 
      patchwork::plot_annotation("$meta.tumour_sample")

    pl_qc = patchwork::wrap_plots(
      CNAqc::plot_peaks_analysis(x, what = 'common', empty_plot = FALSE),
      #CNAqc::plot_peaks_analysis(x, what = 'general', empty_plot = FALSE),
      #CNAqc::plot_peaks_analysis(x, what = 'subclonal', empty_plot = FALSE),
      CNAqc::plot_qc(x),
      CNAqc::plot_CCF(x, assembly_plot = TRUE, empty_plot = FALSE), 
      nrow = 3, 
      heights = c(1,1.5,1)) + 
      patchwork::plot_annotation("$meta.tumour_sample")

    saveRDS(object = x, file = paste0("$prefix", "_qc.rds"))
    saveRDS(object = pl_exp, file = paste0("$prefix", "_plot_data.rds"))
    saveRDS(object = pl_qc, file = paste0("$prefix", "_plot_qc.rds"))

    ggplot2::ggsave(plot = pl_exp, filename = paste0("$prefix", "_data.pdf"), width = 210, height = 297, units="mm", dpi = 200)
    ggplot2::ggsave(plot = pl_qc, filename = paste0("$prefix", "_qc.pdf"), width = 210, height = 297, units="mm", dpi = 200)
    """
}
