process COHORT_SIGNATURES {
    tag "$meta.id"
    label "process_single"
    label "error_retry"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/91/913112a2d9295d35fe97874caf5f500df5a98d6fef1cb5861fd64caa0223a047/data':
        'community.wave.seqera.io/library/r-cnaqc_r-cli_r-dplyr_r-readr_pruned:0fc82bfd06afe6dc' }"

    input:
    tuple val(meta), path(all_fit)
    //tuple val(meta), path(results, stageAs: 'results/*', arity: '0..*')
    //tuple val(meta), path(sparsesignature_assigned, arity: '0..*')
    
    output:
    tuple val(meta), path("*.pdf"), emit: report_cohort_signatures
    tuple val(meta), path("*.rds"), emit: rds_cohort_signatures
    path "versions.yml",            emit: versions


    script:
    def args   = task.ext.args   ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    // Check which inputs are available (they'll be lists)
    //def has_results = results instanceof List ? results.size() > 0 : results as boolean
    //def has_sparse = sparsesignature_assigned instanceof List ? sparsesignature_assigned.size() > 0 : sparsesignature_assigned as boolean

    """
    #!/usr/bin/env Rscript
    print("$all_fit")
    
    library(dplyr)
    library(ggplot2)
    library(ggpubr)
    library(tidyr)
    library(tidyverse)
    # Define color schemes
    sbs_colors = setNames(
      nm = c("SBS1", "SBS17b", "SBS18", "SBS5", "SBS88", "SBS10b", "SBS13", "SBS9", 
             "SBS25", "SBS4", "SBS11", "SBS3", "SBS26", "SBS2", "SBS10a", "SBS31", 
             "Background", "S1", "SBS10c", "SBS10d", "SBS8", "SBS17a", "SBS93", "SBS85"), 
      object = c('#f1696bff', '#8fbd8cff', '#87c7d6ff', '#bac3deff', '#d7bfd9ff', '#a8a2a1ff', 
                 '#cfadb3ff', '#3c609aff', '#9a4564ff', '#fbcb5bff', '#c2b280ff', '#d47e2dff', 
                 '#5f8676ff', 'forestgreen', 'orange', 'brown4', 'grey70', 'darkorchid', 
                 'turquoise4', 'magenta3', 'tan3', 'palegreen1', 'royalblue4', 'olivedrab4')
    )
    
    id_colors = setNames(
      nm = c('ID1', 'ID2', "ID4", "ID5", "ID7", "ID8", "ID9", "ID3", "ID18"), 
      object = c('#0c8281ff', '#f5a55fff', '#7d287eff', '#2e4f4fff', '#c4ddbcff', 
                 '#996869ff', '#daa627ff', 'sienna3', '#bc8f8fff')
    )
    
    dbs_colors <- setNames(
      nm = c("DBS1","DBS2","DBS3","DBS4","DBS5","DBS6","DBS7","DBS8","DBS9","DBS10","DBS11"),
      object = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", 
                 "#a6761d", "#1f78b4", "#b2df8a", "#fb9a99", "#cab2d6")
    )
    
    # Theme function
    my_ggplot_theme <- function () {
      theme_light(base_size=10) +
        theme(legend.key.size=unit(0.3, "cm"),
              panel.background=element_rect(fill="white"),
              axis.text.x=element_text(size=8),
              axis.text.y=element_text(size=8),
              axis.title=element_text(size=10),
              legend.text=element_text(size=8),
              legend.title=element_text(size=10),
              text=element_text(size=10),
              plot.title=element_text(size=12))
    }
    
    # Helper function for exposure extraction
    get_tool_exposure <- function(fit_data, tool, sample_ids=NA) {
      if (tool == "sigprofiler") {
        if (is.null(fit_data)) return(NULL)
        
        fit_data <- as.data.frame(fit_data)
        if (!"Samples" %in% colnames(fit_data)) return(NULL)
        
        fit_data <- fit_data %>%
          tibble::column_to_rownames("Samples") %>%
          mutate(across(everything(), as.numeric))
        
        fit_data_norm <- t(apply(fit_data, 1, function(x) if (sum(x) == 0) x else x / sum(x)))
        if (is.vector(fit_data_norm)) {
          fit_data_norm <- matrix(fit_data_norm, nrow = 1, 
                                  dimnames = list(rownames(fit_data), colnames(fit_data)))
        }
        fit_data <- as.data.frame(fit_data_norm)
      } else if (tool == "sparsesignatures") {
        fit_data <- apply(fit_data, 2, as.numeric)
        n_samples <- nrow(fit_data)
        rownames(fit_data) <- sample_ids
        fit_data <- t(apply(fit_data, 1, function(x) if (sum(x) == 0) x else x / sum(x)))
        fit_data <- as.data.frame(fit_data)
      }
      return(fit_data)
    }
    results_list <- strsplit("$all_fit", " ")[[1]]
    
    # Initialize list for available signatures
    signature_types <- c("SBS96", "ID83")
    availble_sig <- list()
    
    
    check_sigprofiler <- grep(pattern = "SBS96",x = results_list)
    check_sparsesig <- grep(pattern = "cosmic_assigned",x = results_list)
    # has_sparsesig <- grep(pattern = "",x=results_list)
    if (length(check_sigprofiler)!=0){
      has_sigprofiler<-T
    } else{
      has_sigprofiler<-F
    }
    
    if (length(check_sparsesig)!=0){
      has_sparsesig<-T
    } else{
      has_sparsesig<-F
    }
    
    
    # Process SigProfiler results if available
    if (has_sigprofiler) {
      message("Processing SigProfiler results...")
      for (s in signature_types) {
        file_sig <- paste0(s, "/", s, "/Suggested_Solution/COSMIC_", s, 
                           "_Decomposed_Solution/Activities/COSMIC_", s, "_Activities.txt")
        if (file.exists(file_sig)) {
          message("  Found ", s, " activities file")
          sigprofiler_activities <- read.table(file_sig, header = TRUE, sep = "\t")
          sigprofiler_exposures <- get_tool_exposure(fit_data = sigprofiler_activities, 
                                                     tool = "sigprofiler")
          if (!is.null(sigprofiler_exposures)) {
            availble_sig[[s]] <- sigprofiler_exposures %>% 
              tibble::rownames_to_column("sample") %>% 
              pivot_longer(
                cols = !starts_with("sample"),
                names_to = "signature",
                values_to = "value"
              ) %>% 
              mutate(context = s, tool = "SigProfiler")
          }
        } else {
          message("  File for context ", s, " does not exist, skipping...")
        }
      }
      d_sigprofiler_all_long <- do.call("bind_rows", availble_sig)
      # Create cohort-wide exposure plot
      all_cohort_exp <- d_sigprofiler_all_long %>% 
        ggplot(aes(x = sample, y = value, fill = signature)) +
        geom_col() +
        labs(
          x = "Sample ID",
          y = "Exposure",
          fill = "Signature"
        ) +
        my_ggplot_theme() +
        scale_fill_manual(values = c(sbs_colors, id_colors, dbs_colors)) +
        xlab("") +
        theme(panel.spacing.x = unit(0.2, "lines"),
              legend.box = "horizontal",
              panel.grid.minor = element_blank()) +
        facet_wrap(~context, nrow = 3, strip.position = "right") +
        guides(col = guide_legend(nrow = 3)) +
        ggtitle(label = paste0("${meta.dataset}"), 
                subtitle = "Identified signatures in cohort per context")
      
      # Summary of total signatures per sample
      summary_tot_sig <- d_sigprofiler_all_long %>% 
        filter(value != 0) %>% 
        group_by(sample, context) %>% 
        summarise(n = n(), .groups = "drop") %>% 
        ggplot(aes(x = context, y = n)) +
        geom_boxplot() +
        my_ggplot_theme() +
        labs(
          y = "# of signatures",
          x = "Context"
        ) +
        ggtitle("Number of identified signatures per sample")
      
      # Exposure class analysis
      median_value_df <- d_sigprofiler_all_long %>% 
        filter(value != 0) %>% 
        group_by(sample, context) %>% 
        summarise(n = n(), .groups = "drop") %>% 
        group_by(context) %>% 
        summarise(median_signatures = median(n))
      
      exp_classes <- d_sigprofiler_all_long %>% 
        left_join(median_value_df, by = "context") %>%
        group_by(signature) %>% 
        mutate(mean_exposure = mean(value)) %>%
        group_by(context) %>% 
        mutate(mean_exposure_class = 1/median_signatures) %>% 
        mutate(exposure_classes = case_when(
          mean_exposure >= mean_exposure_class ~ "high exposure",
          mean_exposure < mean_exposure_class & mean_exposure >= 0.05 ~ "medium exposure",
          TRUE ~ "low exposure"
        )) %>% 
        filter(exposure_classes != "low exposure") %>% 
        mutate(signature = reorder(signature, -mean_exposure)) %>% 
        select(signature, context, mean_exposure, exposure_classes) %>% 
        distinct() %>% 
        ggplot(aes(x = signature, y = mean_exposure, fill = signature)) +
        geom_col() +
        facet_wrap(~context, scales = "free") +
        scale_fill_manual(values = c(sbs_colors, id_colors, dbs_colors)) +
        my_ggplot_theme() +
        labs(
          y = "Mean exposure",
          x = "Signature"
        ) +
        ggtitle("Mean exposure for most prevalent signatures")
      
      # Cohort prevalence statistics
      n_samples_cohort <- length(unique(d_sigprofiler_all_long\$sample))
      cohort_stats1 <- d_sigprofiler_all_long %>% 
        filter(value != 0) %>% 
        group_by(signature) %>% 
        mutate(n_signatures = n_distinct(sample)) %>% 
        mutate(perc_cohort = n_signatures/n_samples_cohort*100) %>% 
        select(context, signature, perc_cohort) %>% 
        unique() %>% 
        filter(perc_cohort > 10) %>%
        group_by(context) %>% 
        mutate(signature = reorder(signature, -perc_cohort)) %>% 
        ggplot(aes(x = signature, y = perc_cohort, fill = signature)) +
        geom_col(width = 1, color = "white", position = "dodge") +
        facet_wrap(~context, scales = "free") +
        scale_fill_manual(values = c(sbs_colors, id_colors, dbs_colors)) +
        my_ggplot_theme() +
        labs(
          y = "% in cohort",
          x = "Signature"
        ) +
        ggtitle("Most prevalent signatures in cohort per context")
      
    } else {
      all_cohort_exp <- ggplot()
      summary_tot_sig <- ggplot()
      cohort_stats1 <- ggplot()
      exp_classes <- ggplot()
    }
    
    
    # Process SparseSignature results if available
    if (has_sparsesig) {
      message("Processing SparseSignature results...")
      file_sig <- readRDS(paste0("$prefix","_cosmic_assigned.rds"))
      sim_matrix <-file_sig[["similarity_matrix"]]
      
      df <- as.data.frame(sim_matrix)
      
      # Add row names as a column
      df <- df %>%
        rownames_to_column(var = "Sample")
      
      # Pivot to long format
      df_long <- df %>%
        pivot_longer(
          cols = -Sample,
          names_to = "Signature",
          values_to = "Value"
        ) %>% 
        mutate(assignment_class=case_when(Value>=0.95~"Perfect Match",
                                          Value<0.95 & Value>=0.85~"High Cosine",
                                          Value<0.4~"Low Cosine",
                                          TRUE~"Medium Cosine"))
      
      # Plot with geom_tile
      plot_assignment_sparse <- ggplot(df_long, aes(x = Signature, y = Sample, fill = assignment_class)) +
        geom_tile() +
        scale_fill_manual(values=c("Perfect Match"="forestgreen",
                                   "High Cosine"="darkolivegreen4",
                                   "Medium Cosine"="darkolivegreen3",
                                   "Low Cosine"="gainsboro"))+
        ggtitle(label = "Assignment of SparseSignatures signatures to COSMIC Catalogue")+
        my_ggplot_theme()+
        theme(
          axis.text.x = element_text(angle = 90, hjust = 1)
        ) +
        labs(
          x = "SBS Signature",
          y = "DeNovo Signature",
          fill = "Cosine Similarity"
        )
    } else {
      plot_assignment_sparse <- ggplot()
    }
    
    # Combine all signature data
    
    
    
    # Combine plots
    top <- ggarrange(
      all_cohort_exp,
      ncol = 1
    )
    
    middle <- ggarrange(
      cohort_stats1,
      summary_tot_sig,
      ncol = 2,
      widths = c(2, 1),
      legend = "none"
    )
    
    bottom <- ggarrange(
      exp_classes,
      plot_assignment_sparse,
      ncol = 1,
      legend = "none"
    )
    
    final_report <- ggarrange(
      top,
      middle,
      bottom,
      ncol = 1,
      heights = c(2, 1, 3),
      common.legend = TRUE,
      legend = "none"
    )
    
    # Save outputs
    saveRDS(final_report, file = paste0("$prefix", "_cohort_signature_report.rds"))
    ggplot2::ggsave(plot = final_report, filename = paste0("$prefix", "_cohort_signature_report.pdf"), 
                    height = 210, width = 210, units = "mm", dpi = 200)
    # Export versions
    f <- file("versions.yml", "w")
    dplyr_version <- sessionInfo()\$otherPkgs\$dplyr\$Version
    tidyr_version <- sessionInfo()\$otherPkgs\$tidyr\$Version
    ggplot2_version <- sessionInfo()\$otherPkgs\$ggplot2\$Version
    ggpubr_version <- sessionInfo()\$otherPkgs\$ggpubr\$Version
    writeLines(paste0('"', "$task.process", '"', ":"), f)
    writeLines(paste("    dplyr:", dplyr_version), f)
    writeLines(paste("    tidyr:", tidyr_version), f)
    writeLines(paste("    ggplot2:", ggplot2_version), f)
    writeLines(paste("    ggpubr:", ggpubr_version), f)
    close(f)    
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_cohort_signature_report.pdf
    touch ${prefix}_cohort_signature_report.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dplyr: \$(Rscript -e "cat(as.character(packageVersion('dplyr')))")
        tidyr: \$(Rscript -e "cat(as.character(packageVersion('tidyr')))")
    END_VERSIONS
    """
}