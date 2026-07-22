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
    dplyr::filter(VAF > min_VAF) %>% 
    dplyr::filter(blacklisted == FALSE)
  
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
                               ) {
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

# new plotting functions

plot_peaks_analysis = function(x,
                               empty_plot = TRUE,
                               assembly_plot = TRUE,
                               what = "simple") {
  stopifnot(inherits(x, "cnaqc"))
  
  if (what %in% c('simple', 'common'))
  {
    with_peaks = all(!is.null(x\$peaks_analysis))
    if (!with_peaks) {
      warning("Input does not have peaks, see ?peaks_analysis to run peaks analysis.")
      return(CNAqc:::eplot())
    }
    
    karyotypes = x\$peaks_analysis\$fits %>% names
    
    order_karyotypes = c('1:0', '1:1', '2:0', '2:1', '2:2')

    karyotypes = order_karyotypes
    
    # Plot each one of the fits
    plots = lapply(karyotypes, function(k) {
      if (all(is.null(x\$peaks_analysis\$fits[[k]]\$matching)))
      {
        if (empty_plot)
          return(CNAqc:::eplot())
        else
          return(NULL)
      }
      
      return(suppressWarnings(suppressMessages(plot_peaks_fit(x, k))))
    })
    
    plots = plots[!sapply(plots, is.null)]
    if (length(plots) == 0) {
      cli::cli_alert_warning("Nothing to plot")
      return(CNAqc:::eplot())
    }
    
    # Overall QC
    qc = ifelse(x\$peaks_analysis\$QC == 'PASS', 'forestgreen', 'indianred3')
    
    # Plots assembly
    if (assembly_plot)
      plots = suppressWarnings(suppressMessages(
        ggpubr::ggarrange(
          plotlist = plots,
          nrow = 1,
          ncol = length(plots)
        ) +
          ggplot2::theme(
            plot.title = ggplot2::element_text(color = qc),
            panel.border = ggplot2::element_rect(colour = qc,
                                                 fill = NA)
          )
      ))
    
    return(plots)
  }
  
  if (what %in% c('complex', 'general'))
  {
    with_peaks = all(!is.null(x\$peaks_analysis\$general))
    
    if (!with_peaks) {
      warning("Input does not have peaks, see ?peaks_analysis to run peaks analysis.")
      return(eplot())
    }
    
    return(plot_peaks_fit_general(x))
  }
  
  if (what == 'subclonal')
  {
    with_peaks = all(!is.null(x\$peaks_analysis\$subclonal))
    
    if (!with_peaks) {
      warning("Input does not have peaks, see ?peaks_analysis to run peaks analysis.")
      return(CNAqc:::eplot())
    }
    
    pl = plot_peaks_fit_subclonal(x)
    if (assembly_plot)
      pl = ggpubr::ggarrange(
        plotlist = pl,
        ncol = 1,
        nrow = length(pl),
        common.legend = TRUE,
        legend = 'bottom'
      )
    
    return(pl)
  }
  
}

# Plot a single run results with the standard karyotypes model
plot_peaks_fit = function(x, k) {
  matching = x\$peaks_analysis\$matching_strategy
  
  ranges = x\$peaks_analysis\$matches %>%
    dplyr::filter(karyotype == k) %>%
    dplyr::pull(epsilon)
  
  # Required input values
  mutations = x\$mutations %>%
    dplyr::filter(karyotype == k) %>%
    dplyr::mutate(karyotype = paste0(karyotype, " (", matching, ")")) %>% 
    dplyr::filter(!is.na(QC_PASS))
  
  den = x\$peaks_analysis\$fits[[k]]\$density
  expectation = x\$peaks_analysis\$fits[[k]]\$matching %>%
    dplyr::mutate(karyotype = paste0(karyotype, " (", matching, ")"))
  
  xy_peaks = x\$peaks_analysis\$fits[[k]]\$xy_peaks
  purity_error = x\$peaks_analysis\$purity_error

  karyos = x\$n_karyotype[x\$peaks_analysis\$matches\$karyotype %>% unique]
  weight = x\$n_karyotype[k]/sum(karyos)
  
  # Plots cex for anything that is not the main theme
  cex_opt = getOption('CNAqc_cex', default = 1)
  
  # Add QC info
  QC = x\$peaks_analysis\$matches %>%
    dplyr::filter(karyotype == k) %>%
    dplyr::filter(row_number() == 1) %>%
    dplyr::pull(QC)
  
  qc_color = ifelse(QC == "FAIL", "indianred3", 'forestgreen')
  
  
  title = bquote(bold(.(k)) ~
                   .(paste0(
                     ' (n = ', nrow(mutations), ', ', round(weight * 100, 1),  '%)'
                   )))
  
  # Plot the data
  plot_data =
    ggplot2::ggplot(data = mutations, aes(VAF)) +
    ggplot2::geom_histogram(ggplot2::aes(y = ..density..),
                            binwidth = 0.01,
                            alpha = .3) +
    ggplot2::geom_line(
      data = data.frame(x = den\$x, y = den\$y),
      ggplot2::aes(x = x, y = y),
      size = .3,
      color = 'black'
    ) +
    CNAqc:::my_ggplot_theme() +
    ggplot2::labs(title = title,
                  y = 'KDE',
                  x = "VAF") +
    ggplot2::theme(legend.position = 'bottom')  +
    ggplot2::xlim(-0.01, 1.01) +
    ggplot2::facet_wrap( ~ karyotype) +
    ggplot2::theme(strip.background = ggplot2::element_rect(fill = qc_color))
  
  # Add points for peaks to plot
  plot_data = plot_data +
    ggplot2::geom_point(
      data = xy_peaks,
      ggplot2::aes(
        x = x,
        y = y,
        shape = discarded,
        size = counts_per_bin
      ),
      show.legend = FALSE
    ) +
    ggplot2::scale_shape_manual(values = c(`TRUE` = 1, `FALSE` = 16)) +
    ggplot2::scale_size(range = c(1, 3) * cex_opt)
  
  # Add expectation peaks, and matching colors
  plot_data = plot_data +
    ggplot2::geom_point(
      data = expectation,
      ggplot2::aes(x = x, y = y, color = matched),
      size = 2 * cex_opt,
      shape = 4,
      show.legend = FALSE
    ) +
    ggplot2::annotate(
      geom = 'rect',
      xmin = expectation\$x - expectation\$VAF_tolerance,
      xmax = expectation\$x + expectation\$VAF_tolerance,
      ymin = 0,
      ymax = Inf,
      color = NA,
      alpha = .4,
      fill = 'purple4'
    ) +
    ggplot2::geom_segment(
      data = expectation,
      ggplot2::aes(
        x = x,
        y = y,
        xend = peak,
        yend = y,
        color = matched
      ),
      show.legend = FALSE
    ) +
    ggplot2::annotate(
      geom = 'rect',
      xmin = expectation\$peak - expectation\$epsilon,
      xmax = expectation\$peak + expectation\$epsilon,
      ymin = 0,
      ymax = Inf,
      color = NA,
      alpha = .4,
      fill = 'steelblue'
    ) +
    ggplot2::geom_vline(
      data = expectation,
      ggplot2::aes(xintercept = peak, color = matched),
      size = .7 * cex_opt,
      linetype = 'longdash',
      show.legend = FALSE
    ) +
    ggplot2::scale_color_manual(values = c(`TRUE` = 'forestgreen', `FALSE` = 'red'))
  
  # Annotate the offset number
  plot_data = plot_data +
    ggrepel::geom_text_repel(
      data = expectation %>% filter(!matched),
      ggplot2::aes(
        x = x,
        y = y,
        label = round(offset, 2),
        color = matched
      ),
      nudge_x = 0,
      nudge_y = 0,
      size = 3 * cex_opt,
      show.legend = FALSE
    )
  
  return(plot_data)
  
}

# Plot general peaks analysis
plot_peaks_fit_general = function(x) {
  add_counts = function(w) {
    w\$karyotype = paste0(w\$karyotype, ' (n = ', x\$n_karyotype[w\$karyotype], ')')
    w %>% as_tibble()
  }
  
  analysis = x\$peaks_analysis\$general\$analysis
  n_min = x\$peaks_analysis\$general\$params\$n_min
  epsilon =  x\$peaks_analysis\$general\$params\$epsilon
  expected_peaks = x\$peaks_analysis\$general\$expected_peaks %>% add_counts()
  n_bootstrap = x\$peaks_analysis\$general\$params\$n_bootstrap
  data_peaks = x\$peaks_analysis\$general\$data_peaks %>% add_counts()
  data_densities = x\$peaks_analysis\$general\$data_densities  %>% add_counts()
  
  # plotting
  x\$mutations %>%
    dplyr::filter(!is.na(QC_PASS)) %>% 
    filter(karyotype %in% analysis) %>%
    add_counts() %>%
    ggplot2::ggplot(aes(VAF)) +
    ggplot2::geom_histogram(aes(y = ..density..), binwidth = 0.01, fill = 'gray') +
    ggplot2::facet_wrap( ~ karyotype, scales = 'free_y') +
    CNAqc:::my_ggplot_theme() +
    ggplot2::geom_line(data = data_densities,
                       ggplot2::aes(x = x, y = y),
                       inherit.aes = FALSE) +
    ggplot2::geom_vline(
      data = expected_peaks,
      ggplot2::aes(xintercept = peak, color = matched),
      linetype = 'dashed',
      show.legend = FALSE
    ) +
    ggplot2::geom_point(data = data_peaks,
                        ggplot2::aes(x = x, y = y),
                        inherit.aes = FALSE) +
    ggplot2::scale_color_manual(values = c(`FALSE` = 'indianred3', `TRUE` = 'forestgreen')) +
    ggplot2::scale_fill_manual(values = c(`FALSE` = 'indianred3', `TRUE` = 'forestgreen')) +
    ggplot2::geom_rect(
      data = data.frame(
        xmin = expected_peaks\$peak - epsilon,
        xmax = expected_peaks\$peak + epsilon,
        ymin = 0,
        ymax = Inf,
        matched = expected_peaks\$matched,
        karyotype = expected_peaks\$karyotype
      ),
      ggplot2::aes(
        xmin = xmin,
        xmax = xmax,
        ymin = ymin,
        ymax = ymax,
        fill = matched
      ),
      inherit.aes = FALSE,
      alpha = .3
    ) +
    ggplot2::guides(fill = ggplot2::guide_legend("Matched peak", override.aes = ggplot2::aes(alpha = 1))) +
    ggplot2::labs(
      title = bquote(
        "Generalised peak detection ("
        * n['min'] * ' > ' * .(n_min) * ', ' * epsilon * ' = ' * .(epsilon *
                                                                     100) * '%)'
      ),
      caption = bquote(n['nbootstrap'] * ' = ' * .(n_bootstrap))
    )
}

# Plot subclonal peaks analysis
plot_peaks_fit_subclonal = function(x) {
  expected_peaks = x\$peaks_analysis\$subclonal\$expected_peaks
  data_peaks = x\$peaks_analysis\$subclonal\$data_peaks
  data_densities = x\$peaks_analysis\$subclonal\$data_densities
  decision_table = x\$peaks_analysis\$subclonal\$summary
  n_min = x\$peaks_analysis\$subclonal\$params\$n_min
  n_bootstrap = x\$peaks_analysis\$subclonal\$params\$n_bootstrap
  subclonal_mutations = x\$peaks_analysis\$subclonal\$mutations
  epsilon =  x\$peaks_analysis\$subclonal\$params\$epsilon
  
  plot_model_id = function(segment_id)
  {
    this_model_peaks = expected_peaks %>% filter(segment_id == !!segment_id)
    this_model_ids = this_model_peaks\$model_id %>% unique
    
    rank_models = decision_table %>%
      filter(segment_id == !!segment_id) %>%
      arrange(desc(prop)) %>%
      pull(model_id)
    
    rank_models = c(rank_models, setdiff(this_model_ids, rank_models))
    
    which_best = decision_table %>%
      filter(segment_id == !!segment_id) %>%
      arrange(desc(prop))
    which_best = which_best %>% filter(prop == which_best\$prop[1]) %>% pull(model_id)
    
    strip_colors = rep("gray", rank_models %>% length())
    names(strip_colors) = rank_models
    strip_colors[which_best] = "goldenrod3"
    
    rep_muts = lapply(this_model_ids, function(x) {
      subclonal_mutations %>%
        dplyr::filter(!is.na(QC_PASS)) %>% 
        filter(segment_id == !!segment_id) %>%
        mutate(model_id = x,
               model = ifelse(grepl('->', x), "linear", 'branching'))
    }) %>% Reduce(f = bind_rows)
    
    this_title = decision_table %>%
      filter(segment_id == !!segment_id) %>%
      filter(row_number() == 1) %>%
      select(segment_id, size, clones) %>% unlist() %>% paste(collapse = ' ')
    
    rep_muts %>%
      ggplot2::ggplot()  +
      ggplot2::geom_histogram(aes(x = VAF, y = ..density..),
                              binwidth = 0.01,
                              fill = 'gray') +
      ggplot2::xlim(-0.1, 1.1) +
      CNAqc:::my_ggplot_theme() +
      ggplot2::geom_rect(
        data = this_model_peaks,
        ggplot2::aes(
          xmin = peak - epsilon,
          xmax = peak + epsilon,
          fill = matched
        ),
        inherit.aes = FALSE,
        alpha = .3,
        ymin = 0,
        ymax = Inf
      ) +
      ggplot2::scale_color_manual(values = c(`FALSE` = 'indianred3', `TRUE` = 'forestgreen')) +
      ggplot2::scale_fill_manual(values = c(`FALSE` = 'indianred3', `TRUE` = 'forestgreen')) +
      ggplot2::geom_line(
        data = data_densities %>% filter(segment_id == !!segment_id),
        ggplot2::aes(x = x, y = y),
        inherit.aes = FALSE
      ) +
      ggplot2::geom_point(data = data_peaks %>% filter(segment_id == !!segment_id),
                          ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_vline(
        data = expected_peaks %>% filter(segment_id == !!segment_id),
        ggplot2::aes(
          xintercept = peak,
          linetype = role,
          color = matched
        )
      ) +
      ggplot2::facet_wrap( ~ factor(model_id, levels = rank_models)) +
      ggplot2::labs(title = this_title)
    
  }
  
  decision_table\$segment_id %>% unique() %>% lapply(plot_model_id)
  
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
