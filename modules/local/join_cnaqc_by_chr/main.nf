process JOIN_CNAQC_BY_CHR {
    tag "$meta.id"
    label "process_low"
    label "error_retry"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-cnaqc%3A1.1.2--r44hdfd78af_0':
        'biocontainers/r-cnaqc:1.1.2--r44hdfd78af_0' }"

    input:
    tuple val(meta), path(rds_list, stageAs: 'input*.rds'), val(tumour_samples)

    output:
    tuple val(meta), path("*ALL.rds"), val(tumour_samples),  emit: rds_all,  optional: true
    tuple val(meta), path("*PASS.rds"), val(tumour_samples), emit: rds_pass, optional: true
    path "versions.yml",                                     emit: versions

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "$meta.id"
    def qc_filter = args.qc_filter != null ? args.qc_filter : false
    def keep_original = args!="" && args.keep_original ? "$args.keep_original" : ""

    """
    #!/usr/bin/env Rscript

    library(dplyr)
    library(CNAqc)
    library(tibble)

    # defying new cnaqc functions
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

    prepare_input_data = function(mutations, cna, tumour_purity, indels_lenght) {
      # Input data types and fortification
      stopifnot(is_tibble(mutations) | is.data.frame(mutations))
      stopifnot(is_tibble(cna) | is.data.frame(cna))
      stopifnot(tumour_purity > 0 |
                  tumour_purity <= 1 | !is.na(tumour_purity))

      # Input mutations
      ref_nucleotides = c("A", "C", "T", "G")
      alt_nucleotides = c("A", "C", "T", "G")

      mutations = CNAqc:::fortify_mutation_calls(mutations) %>%
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
        }

      }

      return(list(mutations = mutations, cna_clonal = cna_clonal, cna_subclonal = cna_subclonal, tab = tab))
    }

    check_custom_reference = function(x) {
      required = CNAqc::chr_coordinates_hg19 %>% colnames()

      if(!all(colnames(x) %in% required))
      {
        cli::boxx("Problems with your custom reference", background_col = 'red', col = 'white') %>% cat('\n')

        cli::cli_alert_danger("This type of dataframe should be used, but you are missing columns")
        CNAqc::chr_coordinates_hg19 %>% print()

        cli::cli_abort("Cannot use your custom refernece")
      }
    }

    multisample_init <- function(cnaqc_objs, 
                                 # cna_type = "clonal",
                                 QC_filter = TRUE, 
                                 keep_original = TRUE, 
                                 discard_private = FALSE) {
    
      cli::cli_h1("mCNAqcqc - Defining common segments")
      cat('\n')
    
      # perform some checks on the input
      # - it is a list
      # - it contains all CNAqc objects
      # - it is a named list
      # - CNAqc objects all have the 'sample' field
    
      CNAqc:::checking_input(cnaqc_objs)
    
      len = length(cnaqc_objs)
      # cli::cli_alert_info("Selected CNA type: {.field {cna_type}}")
      cli::cli_alert_info("Found {.field {len}} CNAqc objects:")
      cli::cli_ul(names(cnaqc_objs)) 
    
      # create the m_CNAqc object 
    
      cli::cli_h2("Building a {.cls mCNAqcqc} object")
      cat("\n")
    
      # retrive information on breakpoints and define new segments (see the function for better explanation)
      # for the desidered type of mutations
    
      cli::cli_rule("Selecting new segments")
    
      multi_cna <- join_segments(cnaqc_objs = cnaqc_objs, 
                                 # cna_type, 
                                 QC_filter, 
                                 keep_original)
    
      # map the original mutations on the new segments for each segment
    
      cli::cli_rule("Mapping mutations on new segments")
      cat("\n")
    
      multi_mutations <- lapply(names(multi_cna), function(x) {
          set_elements(cnaqc_obj = cnaqc_objs[[x]], 
                     new_cna_list = multi_cna[[x]], 
                     # cna_type = cna_type, 
                     QC_filter = QC_filter#, 
                     # keep_original
                     )
      })
    
      names(multi_mutations) = names(multi_cna)
    
      if(discard_private == TRUE) {
    
        # take only mutations on identical positions across all the samples
        cat("\n")
        cli::cli_rule("Collecting only mutations on shared positions")

        shared_mut = lapply(multi_mutations, function(x) {
          x %>%
            dplyr::mutate(pos = paste(chr, from, to, sep = ":"))
        }) %>% dplyr::bind_rows(.)

        n_samples = shared_mut\$Indiv %>% unique() %>% length()
        tot_n_mut = nrow(shared_mut)
        shared_mut = shared_mut %>%
          dplyr::group_by(pos) %>%
          dplyr::filter(n() == n_samples) %>%
          dplyr::ungroup() %>%
          dplyr::mutate(pos = NULL)

        final_n_mut = nrow(shared_mut)
        removed_n_mut = tot_n_mut - final_n_mut

        cli::cli_alert_info(
          c(
            "Found {.val {tot_n_mut}} mutations mapping common segments across {.val {n_samples}} samples. \n",
            "Removing {.val {removed_n_mut}} mutations on not shared positions, keeping {.val {final_n_mut}} mutations"
          )
        )
        cat("\n")

        multi_mutations = lapply(names(multi_mutations), function(x) {
          shared_mut %>%
            filter(Indiv == x)
        })
        names(multi_mutations) = names(multi_cna)
      }
    
      # create a list with new cnaqc objects with the new segmentation and the mutations mapped on them
    
      cli::cli_h1("Defining the {.cls mCNAqcqc} object")
      cat("\n")
    
      # creating the output 
      m_cnaqc_res <- list()
    
      # define the output as a m_cnaqc object
      class(m_cnaqc_res) <- "m_cnaqc"
    
      # creating the elements of the mCNAqc object
        # list of CNAqc obj with the new segments
    
      new_segmentation_cnaqc <- lapply(names(multi_cna), function(x) {
        #print(x)
        init(
          mutations = multi_mutations[[x]],
          cna = multi_cna[[x]]\$shared,
          purity = cnaqc_objs[[x]]\$purity,
          sample = x, 
          indels_lenght = cnaqc_objs[[x]]\$indels_lenght
        )})
    
      names(new_segmentation_cnaqc) <- sapply(new_segmentation_cnaqc, function(x) {x\$sample})
    
      # multi_input = lapply(names(multi_cna), function(x) {
      #   #print(x)
      #   shared = init(
      #     mutations = multi_mutations[[x]],
      #     cna = multi_cna[[x]]\$shared,
      #     purity = cnaqc_objs[[x]]\$purity,
      #     sample = x
      #   )
    
      m_cnaqc_res\$cnaqc_obj_new_segmentation <- new_segmentation_cnaqc
    
      if (keep_original == TRUE) {
        # private = multi_mutations[[x]]\$private
        # list(mutations_on_shared = shared,
             # mutations_on_private = private,
             # original = cnaqc_objs[[x]]) %>% return()

        m_cnaqc_res\$original_cnaqc_objc <- cnaqc_objs

      } else {

        original_elements <- lapply(cnaqc_objs, function(x) {names(x)}) %>% unlist() %>% unique()
        conserved_elements <- lapply(m_cnaqc_res\$cnaqc_obj_new_segmentation, function(x) {names(x)}) %>% unlist() %>% unique()

        to_include = setdiff(original_elements, conserved_elements)

        other_info = lapply(cnaqc_objs, function(x) {x[to_include]})
        m_cnaqc_res\$original_additional_info <- other_info

        # other_info = lapply(original, function(o) {
        #   cnaqc_objs[[x]][[o]]
        # })
        # names(other_info) = original
        # list(mutations_on_shared = shared,
        #      # mutations_on_private = private,
        #      original_additional_info = other_info) %>% return()
      }
    
    
      # names(multi_input) <- lapply(multi_input, function(x) {x\$mutations_on_shared\$sample}) %>% unlist()
    
      cli::cli_h2("Creating mCNAqc stats")
    
      m_cnaqc_res = generate_mcnaqc_stats(m_cnaqc_res, cnaqc_objs)
    
      if(exists("m_cnaqc_res", inherits = F) & length(m_cnaqc_res\$cnaqc_obj_new_segmentation) == length(cnaqc_objs)) {
        cli::cli_h1("Ended")
        cli::cli_alert_success("{.cls mCNAqc} object created including all samples")
      }
    
      return(m_cnaqc_res) 
      # the output is a m_cnaqc class object, in which each element of the list is a cnaqc object (one per sample) 
      # with all the classical attributes, but on which it has been performed the joint segmentation and mutations 
      # have been therefore remapped
    }  

    get_segment_info = function(data, chr, sample, new_from, new_to, keep_columns){
      data %>%
        dplyr::filter(sample_id == sample,
                      chr == !!chr,
                      from <= new_from,
                      to >= new_to) %>%
        select(all_of(keep_columns))
    }

    # segment definition

    join_segments = function(cnaqc_objs, 
                             # cna_type, 
                             QC_filter, 
                             keep_original){

      # Row binded segments table (with sample specification)
      x = lapply(cnaqc_objs %>% names(), function(x){

        CNA(cnaqc_objs[[x]]#, 
            # type = cna_type
            ) %>% 
          dplyr::mutate(sample_id = x)

      }) %>%
        do.call(bind_rows, .) %>%
        dplyr::select(sample_id, dplyr::everything())
    
      if(QC_filter == TRUE) {
        x = x %>% 
          filter(QC_PASS == TRUE)
        qc_filter_samples = x\$sample_id %>% unique

        if (length(qc_filter_samples) != length(cnaqc_objs %>% names()))  {
          cli::cli_alert_warning(paste(length(setdiff(names(cnaqc_objs), qc_filter_samples)), 
                               'sample(s) did not have any segment passing the QC, excluding it from the creation of the {.cls mCNAqc} object'))
          cat("\n")
        }
      } 
    
    
      out = lapply(x\$chr %>% unique(), function(chr) {

        cli::cli_alert_info("Iterating on {.field {chr}}")
        cli::cli_h2("Iterating on {.field {chr}}")

        old_segments = sapply(x\$sample_id %>% unique(), function(s) {
          x %>%
            dplyr::filter(sample_id == s) %>%
            dplyr::filter(chr == !!chr) %>%
            dplyr::pull(segment_id) %>%
            unique() %>%
            length()
        })

        cli::cli_alert_info("Number of original segments in individual {.cls CNAqc} objects in {.field {chr}}:")
        cli::cli_ul(paste(names(old_segments), old_segments, sep = " = "))
        cat("\n")

        # Chromosome-specific new breakpoints
        new_breakpoints = c(
          x %>%
            dplyr::filter(chr == !!chr) %>%
            dplyr::pull(from),
          x %>%
            dplyr::filter(chr == !!chr) %>%
            dplyr::pull(to)) %>% 
          unique() %>%
          sort()

        cli::cli_alert_info("Found {.field {length(new_breakpoints)}} breakpoints")        

        #  Separate new breakpoints into segment from and to values
        new_from = new_breakpoints[ !new_breakpoints == dplyr::last(new_breakpoints)] # last element will not be included in the from column 
        new_to = new_breakpoints[!new_breakpoints == dplyr::first(new_breakpoints)] # first element will not be included in the to column 

        # iterate on the new breakpoints to subset the cna piled up 
        lapply(new_breakpoints[-1] %>% seq_along(), function(i) {

          # iterate for each sample 
          lapply(x\$sample_id %>% unique(), function(s) {

            # define which columns must be kept in the new table
            not_wanted = c("from", "to", "length", "size", "segment_id", "chr", "sample_id", "n")
            wanted = setdiff(colnames(x), not_wanted)

            # get the information for the sample in the new segment 
            tmp = get_segment_info(x, 
                                   chr = chr, 
                                   sample = s, 
                                   new_from = new_from[i], 
                                   new_to = new_to[i], 
                                   keep_columns = wanted)

            # do some checking on the result
            if (nrow(tmp) == 0) { # there is no information on copy number on the new segment: insert NA as value of all the columns, except from, to, segment_id and sample_id

              tmp_v2 = rep(NA, ncol(tmp))  
              names(tmp_v2) = colnames(tmp)

              tmp = tmp_v2 %>% 
                tibble::as_tibble_row()
            }

            # create a tibble with the information on the new breakpoints and include the previously retrieved information 
            tidyr::tibble(chr = chr, 
                   from = new_from[i], 
                   to = new_to[i], 
                   sample_id = s) %>% 
              dplyr::bind_cols(tmp) %>% 
              dplyr::mutate(segment_id = paste(chr, from, to, sep = ":")) 

          }) %>% do.call(bind_rows, .)
        }) %>% do.call(bind_rows, .)
      })  %>% do.call(bind_rows, .) # create a unique big tibble with the new segmentation
    
    
    # remove all the segments that are not correctly shared across samples
      remove_segments = out %>% 
        dplyr::filter(is.na(Major)) %>% 
        dplyr::filter(is.na(minor)) %>% 
        dplyr::pull(segment_id) %>% 
        unique()
    
      # if(length(remove_segments) != 0) {

      all_segments = out %>%
        pull(segment_id) %>%
        unique()
    
      keep_segments = setdiff(all_segments, remove_segments)
    
      out_shared = out %>%
        dplyr::filter(segment_id %in% keep_segments)

      cli::cli_alert_info("Shared segments across samples: {.val {out_shared %>% pull(segment_id) %>% unique() %>%  length()}}")
    
      # split the shared cna table by sample id
      out_by_sample = lapply(out_shared\$sample_id %>% unique(), function(x) {
        list(shared = out_shared %>%
          dplyr::filter(sample_id == x))
      })
      names(out_by_sample) = lapply(out_by_sample, function(x) {lapply(x, function(y) {y\$sample_id %>% unique}) %>% unlist()}) %>% unlist() %>% unname()
    
      # select mutations on private segments
    
      if (keep_original == TRUE) {

        cli::cli_alert_warning("param {.arg keep_original} is set to {.var TRUE}. Original CNAqc object will be kept.")

        # out_private = out %>%
        #   dplyr::filter(segment_id %in% remove_segments) %>%
        #   dplyr::filter(!is.na(minor) & !is.na(Major))

        # split the private cna table by sample id

        # out_private_by_sample = lapply(out_private\$sample_id %>% unique(), function(x) {
        #   if (nrow(out_private) == 0) {
        #     out_private %>%
        #       dplyr::add_row(chr = NA)
        #   } else {
        #     out_private %>%
        #       dplyr::filter(sample_id == x)
        #   }
        # })
        # 
        # if (length(out_private_by_sample) == 0) {
        #   lapply(x\$sample_id %>% unique(), function(l) {
        #     out_private_by_sample[[l]] = rep(NA, ncol(out_private))
        #   })
        # } else {
        #   names(out_private_by_sample) = lapply(out_private_by_sample, function(x) {
        #     x\$sample_id %>% unique
        #   }) %>% unlist()
        # }
        # 
        # out_by_sample = lapply(names(out_by_sample), function(x) {
        #   out_by_sample[[x]] = list(shared = out_by_sample[[x]]\$shared, private = out_private_by_sample[[x]])
        # })

      } else {
        cat("\n")
        cli::cli_alert_warning(c(
          "param {.arg keep_original} is set to {.var FALSE}.",
          "Found {.val {length(remove_segments)}} not shared segments. Original CNAqc object will not be saved"
          )
        )
      }
    
      names(out_by_sample) = lapply(out_by_sample, function(x) {x\$shared\$sample_id %>% unique()}) %>% unlist()
    
      cli::cli_alert_success("Obtained updated CNA table with shared segments across samples")  
      cat("\n")
    
      return(out_by_sample)
    
      # returns a list with the new segmentation for each sample --> new cna of the m_cnaqc object
    
    }

    # map mutations on new segments
    ## uses CNAqc function "prepare_input_data" to map mutations on the newly defined segments.

    set_elements <- function(cnaqc_obj, 
                             new_cna_list, 
                             # cna_type, 
                             QC_filter #, 
                             # keep_original
                             ) {

      cli::cli_rule(
        crayon::bgCyan(crayon::white(cnaqc_obj\$sample))
      )
    
      cat("\n")
    
      if(class(new_cna_list) != "list") {
        cli::cli_abort(c("Provided a {.cls {class(new_cna_list)}} object", 
                         "x" = "Input must be a list with the new segments"))
      }
    
      initial_mutations = CNAqc::Mutations(cnaqc_obj
                                           # cna = cna_type
                                           )
    
      if(QC_filter == TRUE) {
        initial_mutations = initial_mutations %>% 
          filter(QC_PASS == TRUE)
      }
    
      cli::cli_alert_info("Found {.val {nrow(initial_mutations)}} mutations in the original {.cls CNAqc} object")
      cat("\n")
    
      cli::cli_rule("Mapping mutations on shared segments")
      mutations_shared_segments = prepare_input_data(mutations = initial_mutations, cna = new_cna_list\$shared, tumour_purity = cnaqc_obj\$purity, indels_lenght = cnaqc_obj\$indels_lenght)
      remapped_mut = mutations_shared_segments\$mutations
    
      return(remapped_mut)
    }

    generate_mcnaqc_stats = function(m_cnaqc_res, cnaqc_objs) {
    
      # extract statics 
      n_new_mut <- sapply(m_cnaqc_res\$cnaqc_obj_new_segmentation, function(x) {
        x\$n_mutations
      })
    
      n_or_mut <- sapply(cnaqc_objs, function(x) {
        x\$n_mutations
      })
    
      n_new_cna <- sapply(m_cnaqc_res\$cnaqc_obj_new_segmentation, function(x) {
        x\$n_cna
      })
    
      n_or_cna <- sapply(cnaqc_objs, function(x) {
        x\$n_cna
      })
    
      if (length(n_new_mut) != length(n_or_mut)) {
        missing_samples = setdiff(names(n_or_mut), names(n_new_mut))
        n_new_mut = c(n_new_mut, setNames(NA, missing_samples))
        n_new_cna = c(n_new_cna, setNames(NA, missing_samples))
        all_samples_used_mcnaqc = FALSE
      } else {
        all_samples_used_mcnaqc = TRUE
      }
    
      stats_mcnaqc <- data.frame(
        n_mutations_original_segmentation = n_or_mut,
        n_mutations_new_segmentation = n_new_mut,
        n_cna_original_segmentation = n_or_cna,
        n_cna_new_segmentation = n_new_cna
      )
    
      m_cnaqc_res\$m_cnaqc_stats <- stats_mcnaqc
      m_cnaqc_res\$all_samples_used_mcnaqc = all_samples_used_mcnaqc
    
      return(m_cnaqc_res)
    }

    samples = substr("$tumour_samples", 2, nchar("$tumour_samples")-1)
    samples = strsplit(samples, ", ")[[1]]

    result = lapply(strsplit("$rds_list", " ")[[1]], FUN = function(file){
                readRDS(file)
            })
    names(result) = samples

    result = lapply(result, function(x) {
        x = lapply(x, function(chr) {
            chr\$mutations = chr\$mutations %>% dplyr::rename(Indiv = sample) %>% dplyr::select(-additional_info)
            return(chr)
        })
        return(x)
    })

    chromosomes = lapply(result, names) %>% unlist %>% unique
    
    out_all = lapply(chromosomes, function(cc) {
        mcnaqc_list = lapply(result, function(df) {
            df[[cc]]
        })  
        
        out_all_by_chr = multisample_init(mcnaqc_list,
                                            QC_filter = FALSE,
                                            keep_original = as.logical("$keep_original"),
                                            discard_private = FALSE)
        return(out_all_by_chr)
    }) 
    names(out_all) = chromosomes

    saveRDS(object = out_all, file = paste0("$prefix", "_by_chr_multi_cnaqc_ALL.rds"))

    if (as.logical("$qc_filter") == TRUE){
    
      # select the chr elements
      out_PASS = lapply(chromosomes, function(cc) {
        mcnaqc_list = lapply(result, function(df) {
          df[[cc]]
        })  
        # open the trycatch per chr
        tryCatch(expr = {
          out_PASS_by_chr = multisample_init(mcnaqc_list,
                                                    QC_filter = TRUE,
                                                    keep_original = as.logical("$keep_original"),
                                                    discard_private = FALSE)
          return(out_PASS_by_chr)
        }, error = function(e) {
          print(e)
          print('Not found common segments with QC PASS karyotype, the multi-CNaqc object will be NULL: re-run the pipeline with --filter false')
          out_PASS <<- NULL
        }
        )
      }) 
      names(out_PASS) = chromosomes
      saveRDS(object = out_PASS, file = paste0("$prefix", "_by_chr_multi_cnaqc_PASS.rds"))
    }

    # version export
    f <- file("versions.yml","w")
    dplyr_version <- sessionInfo()\$otherPkgs\$dplyr\$Version
    cnaqc_version <- sessionInfo()\$otherPkgs\$CNAqc\$Version
    writeLines(paste0('"', "$task.process", '"', ":"), f)
    writeLines(paste("    CNAqc:", cnaqc_version), f)
    writeLines(paste("    dplyr:", dplyr_version), f)
    close(f)

    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_multi_cnaqc_ALL.rds
    touch ${prefix}_multi_cnaqc_PASS.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        CNAqc: \$(Rscript -e "cat(as.character(packageVersion('CNAqc')))")
        dplyr: \$(Rscript -e "cat(as.character(packageVersion('dplyr')))")
    END_VERSIONS
    """
}
