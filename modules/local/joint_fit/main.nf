process JOINT_FIT {
  tag "$meta.id"
  container = 'docker://lvaleriani/cnaqc:dev1'
  // container='file:///fast/cdslab/ebusca00/singularity/cdslab.sif'

  input:
    
    tuple val(meta), path(rds_join), path(mobster_best_fits, stageAs:"?/best_fit.rds"), val(tumour_samples)
    // tuple val(datasetID), val(patientID), val(sampleID), path(joint_table), path(mobster_best_fits, stageAs:"?/best_fit.rds")

  output:
    
    tuple val(meta), path("*.rds"), emit: rds //mCNAqc_filtered.rds
    // tuple val(datasetID), val(patientID), val(sampleID), path("$outDir/mCNAqc_filtered.rds")

  script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    //outDir = "subclonal_deconvolution/$datasetID/$patientID"

    //if (!(mobster_best_fits instanceof String)) {
    //  mobster_best_fits = mobster_best_fits.join(",")
    //}

    """
    #!/usr/bin/env Rscript

    # Sys.setenv("VROOM_CONNECTION_SIZE"=99999999)

    library(dplyr)

    patientID = "$meta.patient"
    fits = strsplit(x="$mobster_best_fits", ",") %>% unlist()  # list of mobster fitnames

    check_mutation_id = function(mutations_table) {
      if (!"mutation_id" %in% colnames(mutations_table))
        mutations_table = mutations_table %>% dplyr::mutate(mutation_id=paste(chr,from,to,alt,ref, sep="_"))
      return(mutations_table)
    }

    tail_muts = lapply(fits, function(fit_name) {
      fit_i = readRDS(fit_name)

      data_i = check_mutation_id(fit_i[["data"]])

      data_i %>% 
        dplyr::filter(cluster=="Tail") %>% 
        dplyr::mutate(sample_id=fit_name) %>% 
        dplyr::select(mutation_id, cluster, sample_id)
    }) %>% dplyr::bind_rows() %>% 
      dplyr::group_by(mutation_id) %>% 
      dplyr::summarise(n_samples=dplyr::n()) 
    # dplyr::filter(n_samples == length(fits)) %>%
    # dplyr::pull(mutation_id)

    # Apply filtering based on the value of params_remove_tail
    if ("$params.remove_tail" == "all") {
      muts_to_discard <- tail_muts %>%
        dplyr::filter(n_samples == length(fits)) %>%
        dplyr::pull(mutation_id)
    } else if ("$params.remove_tail" == "once") {
      muts_to_discard <- tail_muts %>%
        dplyr::filter(n_samples >= 1) %>%
        dplyr::pull(mutation_id)
    } else if ("$params.remove_tail" == "never") {
      muts_to_discard <- tail_muts %>%
        dplyr::filter(n_samples < -1000) %>%
        dplyr::pull(mutation_id)
    } else {
      cli::cli_alert_warning("Invalid value for params_remove_tail")
    }

    # Function to discard mutations from the mCNAqc objects
    discard_mutations_from_mCNAqc = function(obj, mutations_list) {
      sample_names = obj[["original_cnaqc_objc"]] %>% names()
      
      get_mCNAqc_types = function(x) {
        names(x) %>% purrr::discard(function(i) i == "m_cnaqc_stats")
      }
      
      types = get_mCNAqc_types(obj)
      
      for (tid in types) {
        sample_names = names(obj[[tid]])
        for (sample_id in sample_names) {
          obj[[tid]][[sample_id]][["mutations"]] = CNAqc::Mutations(obj[[tid]][[sample_id]]) %>% 
            check_mutation_id() %>%
            # dplyr::mutate(unique_id=paste(chr, from, to, ref, alt, sep="_")) %>% 
            dplyr::filter(!mutation_id %in% mutations_list)
        }
      }
      
      return(obj)
    }


    if ( grepl(".rds\$", tolower("$rds_join")) ) {
      obj = readRDS("$rds_join")
      if (class(obj) == "m_cnaqc") {
        obj_filtered = discard_mutations_from_mCNAqc(obj, muts_to_discard)

        saveRDS(obj_filtered, file = paste0("$prefix","_multi_cnaqc_filtered.rds"))
      } else {
        cli::cli_alert_warning("Object of class {class($rds_join)} not supported. Saving the original object.")
        saveRDS(obj, file = paste0("$prefix","_multi_cnaqc_filtered.rds"))
      }
    }

    """
}
