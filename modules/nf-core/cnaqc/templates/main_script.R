#!/usr/bin/env Rscript

# parse_args = function(x) {
#     x = gsub("\\\\[","",x)
#     x = gsub("\\\\]","",x)
#     # giving errors when we have lists like c(xxx, xxx) since it will separate it
#     # args_list = unlist(strsplit(x, ', ')[[1]])
#     args_list = unlist(strsplit(x, ", (?=[^)]*(?:\\\\(|\$))", perl=TRUE))
#     # args_vals = lapply(args_list, function(x) strsplit(x, split=":")[[1]])
#     args_vals = lapply(args_list, function(x) {
#         x_splt = strsplit(x, split=":")[[1]]
#         c(x_splt[1],  paste(x_splt[2:length(x_splt)], collapse=":"))
#     })

#     # Ensure the option vectors are length 2 (key/ value) to catch empty ones
#     args_vals = lapply(args_vals, function(z){ length(z) = 2; z})

#     parsed_args = structure(lapply(args_vals, function(x) x[2]), names = lapply(args_vals, function(x) x[1]))
#     parsed_args[! is.na(parsed_args)]
# }

# opt = list(
#     prefix = ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix')
# )
# args_opt = parse_args('$task.ext.args')
# for ( ao in names(args_opt)) opt[[ao]] = args_opt[[ao]]
# print(opt)

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
  method = "ENTROPY",
  indels_lenght = NULL
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
    stop("matching_epsilon is deprecated - using purity_error = ",
         purity_error)
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

compute_CCF = function(x,
                       karyotypes = c('1:0', '1:1', '2:0', '2:1', '2:2'),
                       muts_per_karyotype = 25,
                       cutoff_QC_PASS = 0.1,
                       method = 'ENTROPY', 
                       min_VAF = 0) {
  stopifnot(inherits(x, 'cnaqc'))
  stopifnot(method %in% c('ENTROPY', "ROUGH"))

  if(any(x\$n_karyotype <= muts_per_karyotype)) warning("Some karyotypes have fewer than", muts_per_karyotype, 'and will not be analysed.')
  nkaryotypes = x\$n_karyotype[x\$n_karyotype > muts_per_karyotype]

  karyotypes = intersect(karyotypes, nkaryotypes %>% names)
  stopifnot(
    karyotypes %in% c('1:0', '1:1', '2:1', '2:0', '2:2')
  )

  # Compute mutation multiplicity
  x\$CCF_estimates = lapply(
    karyotypes,
    function(k)
    {
      if(k %in% c('1:0', '1:1'))
        return(suppressWarnings(mutmult_single_copy(x, k, min_VAF)))

      if(method == "ENTROPY")
        return(suppressWarnings(mutmult_two_copies_entropy(x, k, min_VAF)))
      else
        return(suppressWarnings(mutmult_two_copies_rough(x, k, min_VAF)))
    })
  names(x\$CCF_estimates) = karyotypes

  # Check if there is any null (errors), and remove it
  null_entries = sapply(x\$CCF_estimates, function(x) all(is.null(x)))
  x\$CCF_estimates = x\$CCF_estimates[!null_entries]

  # On extreme cases where there is NO CCF available, we just return x
  if(length(x\$CCF_estimates) == 0) {
    x\$CCF_estimates = NULL
    return(x)
  }

  # Report some stats
  mutations = lapply(x\$CCF_estimates , function(x) x\$mutations)
  mutations = Reduce(dplyr::bind_rows, mutations)

  # pioDisp(
  #   mutations %>%
  #     dplyr::group_by(karyotype, mutation_multiplicity) %>%
  #     dplyr::summarise(assignments = n()) %>%
  #     dplyr::ungroup()
  # )

  # QC the findings
  N = mutations %>%
    dplyr::group_by(karyotype) %>%
    dplyr::summarise(N = n()) %>%
    dplyr::ungroup()

  NA_N = mutations %>%
    dplyr::group_by(karyotype, mutation_multiplicity) %>%
    dplyr::summarise(Unknown = n()) %>%
    dplyr::filter(is.na(mutation_multiplicity)) %>%
    dplyr::select(-mutation_multiplicity) %>%
    dplyr::ungroup()

  QC_table = N %>%
    dplyr::full_join(NA_N, by = 'karyotype') %>%
    dplyr::mutate(
      Unknown = ifelse(is.na(Unknown), 0, Unknown),
      p_Unkown = Unknown/N,
      QC = ifelse(p_Unkown < cutoff_QC_PASS, "PASS", "FAIL"),
      method = method
    )

  if(any(QC_table\$QC == "FAIL"))
  {
    cat('\n')
    cli::cli_h2("Summary CCF assignments. (>{.field {cutoff_QC_PASS*100}%} NAs: not assignable with confidence)")
    print(QC_table)
  }

  for(k in QC_table\$karyotype)
    x\$CCF_estimates[[k]]\$QC_table = QC_table %>% dplyr::filter(karyotype == !!k)

  x
}

mutmult_single_copy = function(x, karyotype, min_VAF) {
  cli::cli_rule("Computing mutation multiplicity for single-copy karyotype {.field {karyotype}}")

  A = as.numeric(strsplit(karyotype, ':')[[1]][1])
  B = as.numeric(strsplit(karyotype, ':')[[1]][2])

  # Karyotype specific mutations - clonal segments
  cl_seg = x\$cna %>%
    dplyr::filter(CCF == 1) %>%
    dplyr::pull(segment_id)

  mutations_k = x\$mutations %>%
    dplyr::filter(VAF > min_VAF) %>%
    filter(blacklisted == FALSE) %>%
    dplyr::filter(karyotype == !!karyotype, segment_id %in% cl_seg) %>%
    dplyr::mutate(
      mutation_multiplicity = 1,
      CCF = CNAqc:::ccf_adjustment_fun(VAF, B, A, x\$purity, mutation_multiplicity)
    )

  return(list(mutations = mutations_k, params = NULL))
}

mutmult_two_copies_entropy = function(x, karyotype, min_VAF) {
  cli::cli_rule(
    "Computing mutation multiplicity for karyotype {.field {karyotype}} using the entropy method."
  )

  A = as.numeric(strsplit(karyotype, ':')[[1]][1])
  B = as.numeric(strsplit(karyotype, ':')[[1]][2])

  # Karyotype specific mutations - clonal segments
  cl_seg = x\$cna %>%
    dplyr::filter(CCF == 1) %>%
    dplyr::pull(segment_id)

  mutations_k = x\$mutations %>%
    dplyr::filter(VAF > min_VAF) %>%
    filter(blacklisted == FALSE) %>%
    dplyr::filter(karyotype == !!karyotype, segment_id %in% cl_seg)

  # Expected VAF for 1 and 2 copies of the mutation
  #
  # Assumption: the aneuploidy state is immediately reached
  # out of a 1:1 state, and therefore we only care about
  # mutations in 1 copy (pre), and 2 copies (post).
  expectation = CNAqc:::expected_vaf_peak(A, B, x\$purity) %>%
    mutate(label = ifelse(mutation_multiplicity == 1, "One copy", "Two copies"))

  med_coverage = median(mutations_k\$DP, na.rm = TRUE)

  cli::cli_alert_info(
    "Expected Binomial peak(s) for these calls (1 and 2 copies): {.value {expectation\$peak}}"
  )

  # =-=-=-=-=-=-=-=-=-=-=-
  # Entropy-derived heuristic for the detection of points
  # that are difficult to assign
  # =-=-=-=-=-=-=-=-=-=-=-
  # We build 2 template Binomial densities to capture:
  #
  # - Bin(p1, n), events before aneuploidy
  # - Bin(p2, n), events after aneuploidy
  #
  # In both cases we take as overal number of trials (n)
  # the median coverage, and use for the success parameters
  # p1 and p2 the expected peaks as of ASCAT equation.
  #
  # Assumptions:
  # - overdispersion is small to justify a Binomial instead
  #   of a Beta-Binomial model;
  # - trials are well-represented with the median coverage;
  p_1 = expectation\$peak[1]
  p_2 = expectation\$peak[2]

  n = ceiling(med_coverage)

  # Bin(p1, n) and Bin(p2, n)
  d_1 = CNAqc:::binomial_density(p_1, n, N_bins = 1000)
  d_2 = CNAqc:::binomial_density(p_2, n, N_bins = 1000)

  # Then we obtain the Binomial quantile ranges for these
  # two distributions, which we use to consider only assingments
  # that have a minimum probability support
  rg_1 = CNAqc:::binomial_quantile_ranges(p_1, n, quantile_left = 0.01, quantile_right = 0.99)
  rg_2 = CNAqc:::binomial_quantile_ranges(p_2, n, quantile_left = 0.01, quantile_right = 0.99)

  # We want to create a mixture model: pi_1 * Bin(p1, n) + (1 - pi_1) * Bin(p2, n)
  # to model the mixture of those two Binomial distributions. To determine
  # the mixing proportions of this mixture we do some empirical trick of
  # get the number of observations between the two Binomial quantile ranges
  # that we have just computed.
  n_rg_1 = mutations_k %>% filter(VAF > rg_1[1], VAF < rg_1[2]) %>% nrow
  n_rg_2 = mutations_k %>% filter(VAF > rg_2[1], VAF < rg_2[2]) %>% nrow

  # Compute the actual mixing proportions, and re-scale the densities accordingly
  mixing = c(n_rg_1, n_rg_2) / (n_rg_1 + n_rg_2)
  d_1\$mixture_y = d_1\$y * mixing[1]
  d_2\$mixture_y = d_2\$y * mixing[2]

  cli::cli_alert_info("Mixing pre/ post aneuploidy: {.value {round(mixing, 2)}}")

  # Now we need to decide how to assign a point in order to determine the actual mutation
  # multeplicity. We want this to be using the entropy of a 2-class model, and the
  # magnitude of the differential of the entropy
  joint = CNAqc:::entropy_profile_2_class(d_1, d_2)

  # if(any(duplicated(joint))) joint = joint[!duplicated(joint), ]

  # Prifile the entropy via peak detection
  entropy_profile_x = joint\$x
  entropy_profile = joint\$entropy

  input_peakdetection = matrix(cbind(x = entropy_profile_x, y = entropy_profile), ncol = 2)
  colnames(input_peakdetection) = c('x', 'y')

  # Peaks detection with these parameters seems to work often
  peaks =  peakPick::peakpick(
    mat = input_peakdetection,
    neighlim = 1,
    deriv.lim = 0.01,
    peak.min.sd = 0,
    peak.npos = 1
  )

  xy_peaks = input_peakdetection[peaks[, 2], , drop = FALSE] %>%
    as_tibble() %>%
    mutate(x = x,
           y = x)


  if (nrow(xy_peaks) == 0) {
    cli::cli_alert_danger("No peaks detected for CCF computation, will not compute values for this karyotype.")

    return(NULL)
  }

  # Points in the centre where there is a violation of the peaks are the actual points we want
  central = entropy_profile_x[which.max(entropy_profile)]

  # signal = entropy_profile
  #
  # J = joint %>%
  #   dplyr::distinct(x, entropy)
  # entropy_profile_x =
  # entropy_profile = joint\$entropy
  #
  # dy = J\$entropy[-1] - J\$entropy[-length(J\$entropy)]
  # dx = J\$x[-1] - J\$x[-length(J\$x)]
  #
  # dydx = dy/dx
  #
  #   infl <- c(FALSE, diff(diff(dydx)>0)!=0)
  #   points(J\$x[infl ], dydx[infl ], col="blue", pch = 3)
  # abline(v=dydx)
  #
  # mdy = abs(median(dy))
  #
  # lp = rp = which.max(J\$entropy)
  #
  #   repeat{
  #     lp = lp - 1
  #     if(lp == 1 | abs(dy[lp]) > mdy) break
  #   }
  #
  #   repeat{
  #     rp = rp + 1
  #     if(rp == length(J\$entropy) | abs(dy[rp]) > mdy) break
  #   }

  lp = xy_peaks %>% dplyr::filter(x < central) %>% dplyr::arrange(desc(x)) %>% dplyr::filter(row_number() == 1) %>% dplyr::pull(x)
  rp = xy_peaks %>% dplyr::filter(x > central) %>% dplyr::arrange(x) %>% dplyr::filter(row_number() == 1) %>% dplyr::pull(x)

  if (length(lp) == 0 | length(rp) == 0) {
    cli::cli_alert_danger(
      "No suitable range of uncertainty detected for CCF, will not compute values for this karyotype."
    )

    return(NULL)
  }


  cli::cli_alert_info("Not assignamble area: [{.value {lp}}; {.value {rp}}]")

  # Assignemnts based on lp
  mutations_k = mutations_k %>%
    rowwise() %>%
    mutate(
      mutation_multiplicity = case_when(VAF <= lp ~ '1',
                                        VAF > rp ~ '2',
                                        TRUE ~ 'NA'),
      mutation_multiplicity = ifelse(
        !is.na(mutation_multiplicity) & mutation_multiplicity != "NA",
        as.numeric(mutation_multiplicity),
        NA
      ),
      CCF = ifelse(
        !is.na(mutation_multiplicity),
        CNAqc:::ccf_adjustment_fun(VAF, B, A, x\$purity, mutation_multiplicity),
        NA
      )
    ) %>%
    ungroup()

  return(list(
    mutations = mutations_k,
    params = list(
      expectation = expectation,
      joint = joint,
      cuts = c(lp, rp),
      method = 'ENTROPY'
    )
  ))
}

mutmult_two_copies_rough = function(x, karyotype, min_VAF) {
  cli::cli_rule(
    "Computing mutation multiplicity for karyotype {.field {karyotype}} using raw VAF cuts."
  )

  A = as.numeric(strsplit(karyotype, ':')[[1]][1])
  B = as.numeric(strsplit(karyotype, ':')[[1]][2])

  # Karyotype specific mutations - clonal segments
  cl_seg = x\$cna %>%
    dplyr::filter(CCF == 1) %>%
    dplyr::pull(segment_id)

  mutations_k = x\$mutations %>%
    dplyr::filter(VAF > min_VAF) %>%
    filter(blacklisted == FALSE) %>%
    dplyr::filter(karyotype == !!karyotype, segment_id %in% cl_seg)

  # Expected VAF for 1 and 2 copies of the mutation as for the entropy case
  expectation = CNAqc:::expected_vaf_peak(A, B, x\$purity) %>%
    mutate(label = ifelse(mutation_multiplicity == 1, "One copy", "Two copies"))

  med_coverage = median(mutations_k\$DP, na.rm = TRUE)

  cli::cli_alert_info(
    "Expected Binomial peak(s) for these calls (1 and 2 copies): {.value {expectation\$peak}}."
  )

  # =-=-=-=-=-=-=-=-=-=-=-
  # Rough-derived heuristic for the detection of points that are difficult to assign
  # =-=-=-=-=-=-=-=-=-=-=-
  p_1 = expectation\$peak[1]
  p_2 = expectation\$peak[2]

  n = ceiling(med_coverage)

  # We get quantiles as for the entropy
  rg_1 = CNAqc:::binomial_quantile_ranges(p_1, n, quantile_left = 0.01, quantile_right = 0.99)
  rg_2 = CNAqc:::binomial_quantile_ranges(p_2, n, quantile_left = 0.01, quantile_right = 0.99)

  # We create the mixture model: pi_1 * Bin(p1, n) + (1 - pi_1) * Bin(p2, n)
  # as with the entropy
  n_rg_1 = mutations_k %>% filter(VAF > rg_1[1], VAF < rg_1[2]) %>% nrow
  n_rg_2 = mutations_k %>% filter(VAF > rg_2[1], VAF < rg_2[2]) %>% nrow

  # So the algebraic midpoint is (p_2 - p_1)/2, we instead split |p_1-p_2|
  # proportionally to n_rg_1 and n_rg_2, normalised
  mixing = c(n_rg_1, n_rg_2) / (n_rg_1 + n_rg_2)
  t_split = p_1 + (p_2 - p_1) * mixing[1]

  cli::cli_alert_info(
    "Mutations per peak: n = {.value {n_rg_1}}, n = {.value {n_rg_2}}. The hard cut is t = {.value {t_split}}."
  )

  # Assignemnts based on lp
  mutations_k = mutations_k %>%
    rowwise() %>%
    mutate(
      mutation_multiplicity = case_when(VAF <= t_split ~ '1',
                                        VAF > t_split ~ '2',
                                        TRUE ~ 'NA'),
      mutation_multiplicity = as.numeric(mutation_multiplicity),
      CCF = ifelse(
        !is.na(mutation_multiplicity),
        CNAqc:::ccf_adjustment_fun(VAF, B, A, x\$purity, mutation_multiplicity),
        NA
      )
    ) %>%
    ungroup()


  return(list(
    mutations = mutations_k,
    params = list(
      expectation = expectation,
      cuts = t_split,
      method = 'ROUGH'
    )
  ))
}

# Script #####


SNV = readRDS("$snv_rds") %>%
  purrr::pluck("$tumour_sample", "mutations") %>%
  dplyr::mutate(mutation_id = paste(chr,from,to,ref,alt,sep = ':'))

CNA = readRDS("$cna_rds")

x = init(mutations = SNV,
        cna = CNA\$segments,
        purity = CNA\$purity,
        sample = "$tumour_sample",
        ref = opt[["genome"]], 
        indels_lenght = opt[["indels_lenght"]])

x = analyze_peaks(x,
                  matching_strategy = opt[["matching_strategy"]],
                  min_absolute_karyotype_mutations = as.numeric(opt[["min_absolute_karyotype_mutations"]]),
                  purity_error = as.numeric(opt[["purity_error"]])
                  )

x = compute_CCF(x,
                muts_per_karyotype = as.numeric(opt[["muts_per_karyotype"]])
)

# this is needed in order to plot the results without the 0 VAF mutations
tmp_x <- x
mut <- CNAqc::Mutations(tmp_x) %>%
  dplyr::filter(VAF > 0)
tmp_x\$mutations <- mut

pl = ggpubr::ggarrange(
  CNAqc::plot_data_histogram(tmp_x, which = 'VAF'),
  CNAqc::plot_data_histogram(tmp_x, which = 'DP'),
  CNAqc::plot_data_histogram(tmp_x, which = 'NV'),
  CNAqc::plot_data_histogram(tmp_x, which = 'CCF'),
  ncol = 2,
  nrow = 2
)

pl_exp = ggpubr::ggarrange(
  plotlist = list(CNAqc::plot_gw_counts(tmp_x),
                  CNAqc::plot_gw_vaf(tmp_x, N = 10000),
                  CNAqc::plot_gw_depth(tmp_x, N = 10000),
                  CNAqc::plot_segments(tmp_x),
                  pl),
  nrow = 5,
  heights = c(.5,.5,.5,1,5)
)
pl_exp = ggpubr::annotate_figure(pl_exp, top = ggpubr::text_grob("$tumour_sample", size = 14))

pl_qc = ggpubr::ggarrange(
  plotlist = list(
    CNAqc::plot_peaks_analysis(tmp_x, what = 'common', empty_plot = FALSE),
    CNAqc::plot_qc(tmp_x),
    CNAqc::plot_CCF(tmp_x, assembly_plot = TRUE, empty_plot = FALSE)),
  nrow = 3,
  heights = c(1,1.5,1))
pl_qc = ggpubr::annotate_figure(pl_qc, top = ggpubr::text_grob("$tumour_sample", size = 14))

saveRDS(object = x, file = paste0(opt[["prefix"]], "_qc.rds"))
saveRDS(object = pl_exp, file = paste0(opt[["prefix"]], "_data_plot.rds"))
saveRDS(object = pl_qc, file = paste0(opt[["prefix"]], "_qc_plot.rds"))

ggplot2::ggsave(plot = pl_exp, filename = paste0(opt[["prefix"]], "_data.pdf"), width = 210, height = 297, units="mm", dpi = 200)
ggplot2::ggsave(plot = pl_qc, filename = paste0(opt[["prefix"]], "_qc.pdf"), width = 210, height = 297, units="mm", dpi = 200)

# version export
f <- file("versions.yml","w")
dplyr_version <- sessionInfo()\$otherPkgs\$dplyr\$Version
cnaqc_version <- sessionInfo()\$otherPkgs\$CNAqc\$Version
writeLines(paste0('"', "$task.process", '"', ":"), f)
writeLines(paste("    CNAqc:", cnaqc_version), f)
writeLines(paste("    dplyr:", dplyr_version), f)
close(f)
