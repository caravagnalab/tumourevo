#!/usr/bin/env Rscript

# parse arguments
parse_args <- function(x){
  args_list <- unlist(strsplit(x, ' ?--')[[1]])[-1]
  args_vals <- lapply(args_list, function(x) scan(text=x, what='character', quiet = TRUE))

  # Ensure the option vectors are length 2 (key/ value) to catch empty ones
  args_vals <- lapply(args_vals, function(z){ length(z) <- 2; z})

  parsed_args <- structure(lapply(args_vals, function(x) x[2]), names = lapply(args_vals, function(x) x[1]))
  parsed_args[! is.na(parsed_args)]
}

# Set defaults and classes

opt <- list(
  prefix = ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix'),
  genome = 'GRCh38',
  genome_coords = NULL,
  karyotypes = c("1:0", "1:1", "2:0", "2:1", "2:2"),
  min_absolute_karyotype_mutations = 100,
  matching_strategy = 'rightmost',
  min_karyotype_size = 0,
  p_binsize_peaks = 0.005,
  matching_epsilon = NULL,
  purity_error = 0.05,
  VAF_tolerance = 0.015,
  n_bootstrap = 1,
  kernel_adjust = 1,
  KDE = TRUE,
  starting_state_subclonal_evolution = "1:1",
  cluster_subclonal_CCF = FALSE,
  min_VAF = 0,
  muts_per_karyotype = 25,
  cutoff_QC_PASS = 0.1,
  method = "ENTROPY"
)
opt_types <- lapply(opt, class)

# Apply parameter overrides

args_opt <- parse_args('$task.ext.args')
for ( ao in names(args_opt)){
  if (! ao %in% names(opt)){
    stop(paste("Invalid option:", ao))
  }else{

    # Preserve classes from defaults where possible
    if (! is.null(opt[[ao]])){
      args_opt[[ao]] <- as(args_opt[[ao]], opt_types[[ao]])
    }
    opt[[ao]] <- args_opt[[ao]]
  }
}

# load libraries
library(dplyr)
library(CNAqc)
library(tibble)
library(gridExtra)
library(ggplot2)

# ---------------------------------------------------------------------------------
# NEW VERSION FUNCTIONS ---------------------------
# ---------------------------------------------------------------------------------

# new init
init = function(mutations, snvs = NULL, cna, purity, sample = "MySample", ref = "GRCh38", genome_coords = NULL, indels_lenght = NULL) {
  cli::cli_h1("CNAqc - CNA Quality Check")
  cat('\n')

  if(!is.null(snvs))
  {
    if(!is.null(mutations))
    {
      cli::boxx("Parameter `snvs` has been deprecated, cannot use it with `mutations`", col = 'red', margin = 3)

      cli::cli_abort("Avoid using `snvs` if you use `mutations")
    }
    else
    {
      cli::boxx("Parameter `snvs` has been deprecated, using it as `mutations",
                col = 'red',
                margin = 3) %>%
        cat('\n')

      mutations = snvs
    }

  }

  # Output
  fit = list()
  class(fit) <- "cnaqc"

# Sample name
  fit\$sample = sample

  # Reference genome
  fit\$reference_genome = ref
  cli::cli_alert_info("Using reference genome coordinates for: {.field {ref}}.")
  
  if(!is.null(genome_coords) & !ref %in% c("hg19", "GRCh37", "hg38", "GRCh38", "mm10", "GRCm38")) {
    fit\$genomic_coordinates = genome_coords
    cli::cli_alert_info("Using custom genome coordinates for: {.field {ref}}.")
  }
  
  # Parse input
  input = prepare_input_data(mutations, cna, purity, indels_lenght = indels_lenght)

  # Remove CNA segments with NA Major/minor

  fit\$mutations = input\$mutations
  fit\$cna = input\$cna_clonal %>%
    dplyr::left_join(input\$tab, by = 'segment_id') %>% as_tibble()
  fit\$cna_subclonal = input\$cna_subclonal %>% as_tibble()
  fit\$has_subclonal_CNA = !all(is.null(input\$cna_subclonal))

  # Counts data
  if(!fit\$has_subclonal_CNA)
    fit\$n_mutations = nrow(fit\$mutations)
  else
    fit\$n_mutations = nrow(fit\$mutations)  + sapply(fit\$cna_subclonal\$mutations, nrow) %>% sum()

  fit\$n_cna_clonal = nrow(fit\$cna)
  fit\$n_cna_subclonal = ifelse(is.null(fit\$cna_subclonal), 0, nrow(fit\$cna_subclonal))
  fit\$n_cna = fit\$n_cna_clonal + fit\$n_cna_subclonal

  fit\$n_karyotype = sort(table(fit\$mutations\$karyotype), decreasing = T)
  fit\$purity = purity

  # Segments length (clonal)
  genome_segs_length = fit\$cna %>%
    dplyr::group_by(Major, minor) %>%
    dplyr::summarise(L = sum(length), .groups = 'drop') %>%
    dplyr::mutate(karyotype = paste0(Major, ':', minor)) %>%
    dplyr::arrange(desc(L))

  fit\$l_karyotype = genome_segs_length\$L
  names(fit\$l_karyotype ) = genome_segs_length\$karyotype

  tab_ploidy = fit\$cna %>%
    dplyr::group_by(minor, Major) %>%
    dplyr::summarise(n = sum(length)) %>%
    dplyr::arrange(desc(n)) %>%
    dplyr::mutate(karyotype = paste0(Major, ':', minor)) %>%
    dplyr::ungroup()

  fit\$ploidy = as.numeric(tab_ploidy\$minor[1]) + as.numeric(tab_ploidy\$Major[1])
  fit\$most_prevalent_karyotype = paste0(tab_ploidy\$Major[1], ':',  tab_ploidy\$minor[1])
  fit\$basepairs_by_karyotype = tab_ploidy

  fit\$most_mutations_karyotype = names(fit\$n_karyotype)[1]
  
  fit\$indels_lenght = indels_lenght
  
  fit
}

# prepare input data
prepare_input_data = function(mutations, cna, tumour_purity, indels_lenght) {
  # Input data types and fortification
  stopifnot(is_tibble(mutations) | is.data.frame(mutations))
  stopifnot(is_tibble(cna) | is.data.frame(cna))
  stopifnot(tumour_purity > 0 |
              tumour_purity <= 1 | !is.na(tumour_purity))

  # Input mutations
  ref_nucleotides = c("A", "C", "T", "G")
  alt_nucleotides = c("A", "C", "T", "G")

  mutations =  CNAqc:::fortify_mutation_calls(mutations) %>%
    mutate(
      type = ifelse(
        (ref %in% ref_nucleotides) & (alt %in% alt_nucleotides),
        'SNV',
        'indel'
      )
    )
  
  if(!is.null(indels_lenght)) {
    mutations = mutations %>% 
      mutate(blacklisted = case_when(
        type == 'SNV' ~ FALSE, 
        (type == 'indel' &  (nchar(alt) >= indels_lenght | nchar(ref) >= indels_lenght)) ~ TRUE, 
        (type == 'indel' &  (nchar(alt) < indels_lenght | nchar(ref) < indels_lenght)) ~ FALSE))
  } else {
    mutations = mutations %>% 
      mutate(blacklisted = FALSE)
  }

  all_mutations = mutations

  # Check specific mappability features for driver
  driver_mutations = NULL
  if(all(c("is_driver", "driver_label") %in% colnames(all_mutations)))
  {
    driver_mutations = all_mutations %>% filter(is_driver)

    if(nrow(driver_mutations) == 0) driver_mutations = NULL
  }

  if(!is.null(driver_mutations))
  {
    cli::cli_alert_success("Found annotated driver mutations: {.field {driver_mutations\$driver_label}}.")
  }

  nsnvs = (all_mutations\$type == "SNV") %>% sum()
  nindel = nrow(all_mutations) - nsnvs
  psnvs = round(nsnvs/(nsnvs + nindel) * 100)
  cli::cli_alert_success("Fortified calls for {.field {nrow(all_mutations)}} somatic mutations: {.field {nsnvs}} SNVs ({.field {psnvs}%}) and {.field {nindel}} indels.")

  # CNAs - clonal and subclonal
  cna_all = CNAqc:::fortify_CNA_segments(cna)

  cna_clonal = cna_all\$clonal %>% mutate(segment_id = paste(chr, from, to, Major, minor, CCF, sep = ':'))

  cna_subclonal = cna_all\$subclonal
  if(!is.null(cna_subclonal)) cna_subclonal = cna_subclonal %>%
    mutate(segment_id = paste(chr, from, to, Major, minor, CCF, sep = ':'))

  ncnacl = cna_clonal %>% nrow()
  ncnasbcl = ifelse(is.null(cna_subclonal), 0, cna_subclonal %>% nrow())

  cli::cli_alert_success("Fortified CNAs for {.field {ncnacl + ncnasbcl}} segments: {.field {ncnacl}} clonal and {.field {ncnasbcl}} subclonal.")

  # Mapping mutations
  # cli::cli_alert_info("Mapping {.field {nrow(snvs)}} mutations")
  #   paste0("Input ",
  #   'n = ',
  #   nrow(snvs),
  #   " mutations for ",
  #   nrow(cna),
  #     " CNA segments (",
  #     ncnacl,
  #     " clonal, ",
  #     ncnasbcl,
  #     " subclonal)"
  #   )
  # )

  # Mapping mutations to clonal segments
  mutations = CNAqc:::map_mutations_to_clonal_segments(mutations, cna_clonal)
  
  # Notify NA (not mapped)
  not_mappable = sum(is.na(mutations\$karyotype))

  if(not_mappable > 0)
  {
    # cli::cli_alert_danger('{.field {not_mappable}} mutations cannot be mapped to clonal CNAs and will be removed.')
    mutations = mutations %>% dplyr::filter(!is.na(karyotype))
  }

  nmutations = nrow(mutations)

  # Tabular of mapping per segments (count)
  tab = mutations %>%
    group_by(segment_id) %>%
    summarise(n = n()) %>%
    arrange(desc(n))

  # Stats about mappability
  num_mappable = sum(!is.na(mutations\$karyotype))
  perc_mappable = round(num_mappable / nsnvs * 100)

  cli::cli_alert_success(
    "{.field {num_mappable}} mutations mapped to clonal CNAs."
  )

  # Mapping mutations subclonal segments
  if(!is.null(cna_subclonal))
    {
    cna_subclonal = CNAqc:::map_mutations_to_subclonal_segments(mutations = all_mutations, cna_subclonal)

    if(nrow(cna_subclonal) == 0)
    {
      cna_subclonal = NULL
      cli::cli_alert_info("No mutations mapped to subclonal CNAs.")
    }
    else
    {
      nsubcl = sapply(cna_subclonal\$mutations, nrow) %>% sum
      cli::cli_alert_success(
        "{.field {nsubcl}} mutations mapped to subclonal CNAs."
      )
    }
  }

  # Check if we mapped all the drivers
  idfy = function(x) { x %>% dplyr::mutate(id = paste(chr, from, to, ref, alt)) }

  if(!is.null(driver_mutations))
  {
    ids_clonal = mutations %>% idfy %>% pull(id)
    ids_subclonal = NULL

    if(!is.null(cna_subclonal))
      ids_subclonal = cna_subclonal\$mutations %>%
        Reduce(f = dplyr::bind_rows) %>%
        idfy %>% pull(id)

    all_ids = c(ids_clonal, ids_subclonal)
    driver_ids = driver_mutations %>% idfy %>% pull(id)

    missing_drivers = which(!(driver_ids %in% all_ids), arr.ind = TRUE)

    if(length(missing_drivers) > 0)
    {
      missing = driver_mutations\$driver_label[missing_drivers]
      missing = paste("Driver(s): ", paste(missing, collapse = ', '))

      cat("\n")
      cli::boxx(
        "Driver cannot be mapped - out of any segment!",
        header = missing,
        float = "center",
        col = 'white',
        border_col = 'black',
        background_col = 'indianred3'
      ) %>% cat
      cat("\n")

      # cli::boxx(
      #   "Lost driver mutation(s) during mappability (reported below) - out of any segment!",
      #   float = "center",
      #   col = 'white',
      #   background_col = 'red'
      #   ) %>% cat
      #
      # driver_mutations[missing_drivers, ] %>%
      #   dplyr::select(
      #     chr, from, to, ref, alt, NV, DP, VAF,
      #     is_driver, driver_label
      #   ) %>%
      #   print()
    }

  }

  return(list(mutations = mutations, cna_clonal = cna_clonal, cna_subclonal = cna_subclonal, tab = tab))
}

analyze_peaks = function(x,
                        karyotypes = c('1:0', '1:1', '2:0', '2:1', '2:2'),
                        min_karyotype_size = 0,
                        min_absolute_karyotype_mutations = 100,
                        p_binsize_peaks = 0.005,
                        matching_epsilon = NULL,
                        purity_error = 0.05,
                        VAF_tolerance = 0.015,
                        n_bootstrap = 1,
                        kernel_adjust = 1,
                        matching_strategy = "closest",
                        KDE = TRUE,
                        starting_state_subclonal_evolution = "1:1",
                        cluster_subclonal_CCF = FALSE,
                        min_VAF = 0) {
  
  if (!is.null(matching_epsilon)) {
    stop("matching_epsilon is deprecated - using purity_error = ", purity_error)
    matching_epsilon = purity_error
  }
  
  # Check inputs
  stopifnot(inherits(x, "cnaqc"))
  stopifnot(min_karyotype_size >= 0 & min_karyotype_size < 1)
  stopifnot(min_absolute_karyotype_mutations >= 0)
  stopifnot(p_binsize_peaks > 0 & p_binsize_peaks < 1)
  stopifnot(purity_error > 0 & purity_error < 1)
  stopifnot(is.numeric(kernel_adjust))
  stopifnot(matching_strategy %in% c("closest", "rightmost"))
  
  # Common peaks analysis - they must be in the sample
  cli::cli_h1("Peak analysis: simple CNAs")
  cat("\n")
  
  # filter the mutations using min_VAF and recompute the number of muts per karyotype 
  # and other statistics --> in this way the peak analysis is more precise
  
  # save the mutations with VAF < min_VAF -- they might be uselful after the peak analysis!
  
  blacklisted_muts = x\$mutations %>%
    dplyr::filter(VAF <= min_VAF | blacklisted)
  
  x\$mutations = x\$mutations %>%
    dplyr::filter(VAF > min_VAF)
  
  # recompute the number of mutations mapping on each karyotype
  n_karyo = x\$mutations %>% 
    dplyr::group_by(karyotype) %>% 
    dplyr::summarize( n = dplyr::n())
  
  # change the values of the class
  x\$n_karyotype = setNames(nm = n_karyo\$karyotype, object = n_karyo\$n)
  
  x\$most_mutations_karyotype = n_karyo %>% 
    dplyr::slice_max(n) %>% 
    dplyr::pull(karyotype)

  x = x %>% CNAqc:::analyze_peaks_common(
    karyotypes = karyotypes,
    min_karyotype_size = min_karyotype_size,
    min_absolute_karyotype_mutations = min_absolute_karyotype_mutations,
    p_binsize_peaks = p_binsize_peaks,
    purity_error = purity_error,
    VAF_tolerance = VAF_tolerance,
    n_bootstrap = n_bootstrap,
    kernel_adjust = kernel_adjust
    # min_VAF = min_VAF
  )
  
  x\$cna <- dplyr::left_join(x\$cna,CNAqc:::compute_QC_table(x)\$QC_table %>% filter(type == "Peaks") %>% 
                              dplyr::mutate(QC_PASS = dplyr::if_else(QC == "PASS", TRUE, FALSE)) %>%
                              tidyr::separate(karyotype, into = c("Major", "minor"), sep = ":") %>%
                              mutate(Major = as.numeric(Major), minor = as.numeric(minor)) %>% 
                              dplyr::select(Major, minor, QC_PASS))
  x\$mutations <- dplyr::left_join(x\$mutations %>% ungroup(), CNAqc:::compute_QC_table(x)\$QC_table %>%
                                    filter(type == "Peaks") %>% 
                                    dplyr::mutate(QC_PASS = dplyr::if_else(QC == "PASS", TRUE, FALSE)) %>% 
                                    dplyr::select(karyotype, QC_PASS))
  
  # Generalised peak analysis
  cli::cli_h1("Peak analysis: complex CNAs")
  cat("\n")
  
  w = x\$n_karyotype[!(x\$n_karyotype %>% names() %in% karyotypes)]
  w = w[w > min_absolute_karyotype_mutations]
  
  if(length(w) > 0)
  {
    cli::cli_alert_info(
      "Karyotypes {.field {names(w)}} with >{.field {min_absolute_karyotype_mutations}} mutation(s). Using epsilon = {.field {purity_error}}."
    )
    
    x = x %>% CNAqc:::analyze_peaks_general(
      n_min = min_absolute_karyotype_mutations,
      epsilon = purity_error,
      kernel_adjust = kernel_adjust,
      n_bootstrap = n_bootstrap 
      # min_VAF = min_VAF
    )
    
    x\$peaks_analysis\$general\$summary %>%
      print()
    
    x\$cna <- dplyr::left_join(x\$cna, x\$peaks_analysis\$general\$summary %>% 
                                dplyr::mutate(QC_PASS = dplyr::if_else(prop >= 0.5, TRUE, FALSE)) %>%
                                tidyr::separate(karyotype, into = c("Major", "minor"), sep = ":") %>%
                                mutate(Major = as.numeric(Major), minor = as.numeric(minor)) %>% 
                                dplyr::select(Major, minor, QC_PASS))
    x\$mutations <- dplyr::left_join(x\$mutations, x\$peaks_analysis\$general\$summary %>% 
                                      dplyr::mutate(QC_PASS = dplyr::if_else(prop >= 0.5, TRUE, FALSE)) %>% 
                                      dplyr::select(karyotype, QC_PASS))
  }
  else
    cli::cli_alert_info(
      "No karyotypes with >{.field {min_absolute_karyotype_mutations}} mutation(s). "
    )
  
  cli::cli_h1("Peak analysis: subclonal CNAs")
  cat("\n")
  
  # Subclonal CNAs peak analysis
  if(x\$n_cna_subclonal > 0) {
    x = x %>% CNAqc:::analyze_peaks_subclonal(
      n_min = min_absolute_karyotype_mutations,
      epsilon = purity_error,
      kernel_adjust = kernel_adjust,
      n_bootstrap = n_bootstrap,
      starting_state = starting_state_subclonal_evolution,
      cluster_subclonal_CCF = cluster_subclonal_CCF
    )
    
    if(!is.null(x\$peaks_analysis\$subclonal)) {
      x\$peaks_analysis\$subclonal\$summary %>% print()
    } else {
      cli::cli_alert_info("Subclonal CNAs not analysed with the current parameters.")
      }
  } else {
    cli::cli_alert_info("No subclonal CNAs in this sample.")
    }
  
  x\$mutations = bind_rows(x\$mutations, blacklisted_muts)
  
  return(x)
}

split_by_chromosome = function(x,
                               chromosomes = paste0('chr', c(1:22, 'X', 'Y'))
                               )
{
  stopifnot(inherits(x, 'cnaqc'))

  objs = NULL
  nm = NULL

  for(chr in chromosomes)
  {
    clonal_mutations = x\$mutations %>% filter(chr == !!chr)

    if((clonal_mutations %>% nrow()) == 0) next

    cli::cli_h3(chr)
    cat("\n")
    
    subclonal_mutations = NULL

    if(nrow(x\$cna_subclonal) > 0) subclonal_mutations = x\$cna_subclonal %>%
      filter(chr == !!chr) %>% pull(mutations) %>%
      Reduce(f = bind_rows)
    
    cnas = x\$cna %>% filter(chr == !!chr)
    
    
    if(!is.null(x\$genome_coords) & !x\$reference_genome %in% c("hg19", "GRCh37", "hg38", "GRCh38", "mm10", "GRCm38")) {
      cnaqc_obj = init(
        clonal_mutations %>% bind_rows(subclonal_mutations),
        cna = cnas,
        purity = x\$purity,
        ref = x\$reference_genome,
        genome_coords = x\$genomic_coordinates, 
        indels_lenght = x\$indels_lenght
      )
    } else { 
      cnaqc_obj = init(
        clonal_mutations %>% bind_rows(subclonal_mutations),
        cna = cnas,
        purity = x\$purity,
        ref = x\$reference_genome, 
        indels_lenght = x\$indels_lenght
      )}

    objs = append(objs, list(cnaqc_obj))
    nm = c(nm, chr)
  }

  names(objs) = nm

  return(objs)
  
}


# Script #####

x = readRDS('$cnaqc_rds')

x_by_chr = split_by_chromosome(x)
x_by_chr = lapply(x_by_chr, function(dd) {
    analyze_peaks(dd, 
        matching_strategy = opt[["matching_strategy"]],
        min_absolute_karyotype_mutations = as.numeric(opt[["min_absolute_karyotype_mutations"]]),
        purity_error = as.numeric(opt[["purity_error"]])
    )
})

# now add the plots
# this is needed in order to plot the results without the 0 VAF mutations
tmp_x <- x_by_chr
tmp_x = lapply(tmp_x, function(chr) {
    chr\$mutations <- chr\$mutations %>% 
        dplyr::filter(VAF > 0)  

    new_id = paste(chr\$sample, unique(chr\$mutations\$chr), sep = '_')

    chr\$sample = new_id

    return(chr)
})

# vaf_plt = lapply(tmp_x, function(x) {plot_data_histogram(x, which = 'VAF')})
# vap_plt = ggpubr::ggarrange(plotlist = vaf_plt)

# dp_plt = lapply(tmp_x, function(x) {plot_data_histogram(x, which = 'DP')})
# dp_plt = ggpubr::ggarrange(plotlist = dp_plt)

# nv_plt = lapply(tmp_x, function(x) {plot_data_histogram(x, which = 'NV')})
# nv_plt = ggpubr::ggarrange(plotlist = nv_plt)

pa_plt = lapply(tmp_x, function(cc) {
  plot_peaks_analysis(cc, what = 'common', empty_plot = FALSE) +
    ggtitle(cc\$sample)
})
pa_plt <- marrangeGrob(pa_plt, nrow = 4, ncol = 1, top = NULL)

saveRDS(object = x_by_chr, file = paste0(opt[["prefix"]], "_by_chr_qc.rds"))
saveRDS(object = pa_plt, file = paste0(opt[["prefix"]], "__by_chr_qc_plot.rds"))

ggplot2::ggsave(plot = pa_plt, filename = paste0(opt[["prefix"]], "_by_chr_qc.pdf"), width = 210, height = 297, units="mm", dpi = 200)

# version export
f <- file("versions.yml","w")
dplyr_version <- sessionInfo()\$otherPkgs\$dplyr\$Version
cnaqc_version <- sessionInfo()\$otherPkgs\$CNAqc\$Version
gridextra_version <- sessionInfo()\$otherPkgs\$gridExtra\$Version
writeLines(paste0('"', "$task.process", '"', ":"), f)
writeLines(paste("    CNAqc:", cnaqc_version), f)
writeLines(paste("    dplyr:", dplyr_version), f)
writeLines(paste("    gridExtra:", gridextra_version), f)
close(f)
