process COHORT_SIGNATURES {
    tag "$meta.id"
    label "process_low"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/53/53da5a4364ce10f58c67142c1b4203f19eeff99b533c3e3a310d0a3b37fb46e8/data':
        'community.wave.seqera.io/library/r-ggforce_r-ggnewscale_r-ggraph_r-ggrepel_pruned:5408b55f63d978d5' }"

    input:
    tuple val(meta), path(all_fit)
    
    output:
    tuple val(meta), path("*.pdf"), emit: report_cohort_signatures
    tuple val(meta), path("*.rds"), emit: rds_cohort_signatures
    tuple val(meta), path("*.txt"), emit: table_cohort_signatures
    path "versions.yml",            emit: versions


    script:
    def args   = task.ext.args   ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    #!/usr/bin/env Rscript
    library(tidyverse)
    library(ggplot2)
    library(RColorBrewer)
    library(patchwork)
    # Define color schemes
    COSMIC_color_palette = function(seed=55) {
    COSMIC_sbs <- c(
        "SBS1", "SBS2", "SBS3", "SBS4", "SBS5", "SBS6",
        "SBS7a", "SBS7b", "SBS7c", "SBS7d",
        "SBS8", "SBS9",
        "SBS10a", "SBS10b", "SBS10c", "SBS10d",
        "SBS11", "SBS12", "SBS13", "SBS14", "SBS15", "SBS16",
        "SBS17a", "SBS17b",
        "SBS18", "SBS19", "SBS20", "SBS21",
        "SBS22a", "SBS22b",
        "SBS23", "SBS24", "SBS25", "SBS26", "SBS27", "SBS28",
        "SBS29", "SBS30", "SBS31", "SBS32", "SBS33", "SBS34",
        "SBS35", "SBS36", "SBS37", "SBS38", "SBS39",
        "SBS40a", "SBS40b", "SBS40c",
        "SBS41", "SBS42", "SBS43", "SBS44", "SBS45", "SBS46",
        "SBS47", "SBS48", "SBS49", "SBS50",
        "SBS51", "SBS52", "SBS53", "SBS54", "SBS55", "SBS56",
        "SBS57", "SBS58", "SBS59", "SBS60",
        "SBS84", "SBS85", "SBS86", "SBS87", "SBS88", "SBS89",
        "SBS90", "SBS91", "SBS92", "SBS93", "SBS94", "SBS95",
        "SBS96", "SBS97", "SBS98", "SBS99"
    )
    COSMIC_dbs <- c(
        "DBS1", "DBS2", "DBS3", "DBS4", "DBS5",
        "DBS6", "DBS7", "DBS8", "DBS9", "DBS10",
        "DBS11", "DBS12", "DBS13", "DBS14", "DBS15",
        "DBS16", "DBS17", "DBS18", "DBS19", "DBS20"
    )
    COSMIC_indels <- c(
        "ID1", "ID2", "ID3", "ID4", "ID5", "ID6",
        "ID7", "ID8", "ID9", "ID10", "ID11", "ID12",
        "ID13", "ID14", "ID15", "ID16", "ID17", "ID18"
    )
    catalogues = list(COSMIC_dbs,
                        COSMIC_indels, COSMIC_sbs)
    sigs = catalogues %>% unlist() %>% unique()
    set.seed(seed)
    color_palette <- colors(distinct = length(sigs))
    names(color_palette)<-sigs
    return(color_palette)
    }
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
    } else {
    d_sigprofiler_all_long <- data.frame()
    }


    # Process SparseSignature results if available
    if (has_sparsesig) {
    message("Processing SparseSignature results...")
    file_sig <- readRDS(paste0("$prefix","_cosmic_assigned.rds"))
    sprsesignatures_exposures <- file_sig\$remapped_exposures_prop %>% as.data.frame()
    d_sparsesignatures_all_long <- sprsesignatures_exposures %>% 
        tibble::rownames_to_column("sample") %>% 
        pivot_longer(
        cols = !starts_with("sample"),
        names_to = "signature",
        values_to = "value"
        ) %>% 
        mutate(context = "SBS96", tool = "SparseSignatures")
    cosmic_assign_sparse_sig <- file_sig\$similarity_matrix %>% 
        as.data.frame() %>% 
        rownames_to_column("group") %>% 
        pivot_longer(
        cols = -group,
        names_to = "signature",
        values_to = "value"
        )
    cosmic_assign_sparse_sig_plot <- ggplot(cosmic_assign_sparse_sig, aes(
        x = signature,
        y = group,
        fill = value
    )) +
        geom_tile(color = "white") +
        scale_fill_gradient(
        low = "white",
        high = "red",
        name = "Cosine similarity"
        )+
        theme_minimal() +
        theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank()
        ) +
        labs(
        x = "Signature",
        y = NULL
        )+
        ggtitle("COSMIC assignment of SparseSignatures de novo SBS")
    } else {
        d_sparsesignatures_all_long <- data.frame()
        cosmic_assign_sparse_sig_plot <- ggplot()
    }

    # Combine all signature data
    d_all_signatures <- rbind(d_sigprofiler_all_long,d_sparsesignatures_all_long)
    tot_samples <- d_all_signatures %>% pull(sample) %>% unique() %>% length()
    all_cohort_exp <- d_all_signatures %>%
    ggplot(aes(x = sample, y = value, fill = signature)) +
    geom_col() +
    labs(
        x = "Sample ID",
        y = "Exposure",
        fill = "Signature"
    ) +
    my_ggplot_theme() +
    scale_fill_manual(values = COSMIC_color_palette()) +
    xlab("") +
    theme(panel.spacing.x = unit(0.2, "lines"),
            legend.box = "horizontal",
            panel.grid.minor = element_blank(),axis.text.x = element_blank()) +
    facet_grid(tool~context) +
    ggtitle(label = paste0("Identified signatures in ","${meta.dataset}", " cohort per context"))+
    theme(legend.position = "bottom")+
    guides(
        fill = guide_legend(
        nrow = 2,
        byrow = TRUE
        )
    )



    summary_tot_sig <- d_all_signatures %>%
    filter(value != 0) %>%
    group_by(sample, context, tool) %>%
    summarise(n = n(), .groups = "drop") %>%
    ggplot(aes(x = context, y = n)) +
    geom_boxplot() +
    facet_wrap(~tool) +
    my_ggplot_theme() +
    labs(
        y = "# of signatures",
        x = "Context"
    ) +
    ggtitle( "Number of identified signatures per sample")

    # Exposure class analysis
    median_value_df <- d_all_signatures %>%
    filter(value != 0) %>%
    group_by(sample, context, tool) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(context, tool) %>%
    summarise(median_signatures = median(n))


    exp_classes_all <-d_all_signatures %>%
    filter(value != 0) %>%
    ggplot(aes(
        x = signature,
        y = value,
        fill = signature
    )) +
    geom_boxplot(outlier.alpha = 0.4)+
    facet_grid(tool~ context, scales = "free") +
    scale_fill_manual(values = COSMIC_color_palette())+
    my_ggplot_theme()+
    theme(legend.position = "none",
            axis.text.x = element_text(hjust = 1,angle = 45))+
    labs(
        x = "Signature",
        y = "Tool"
    )+
    ggtitle("Exposure distribution of identified signatures")


    frq_cohort_plot <- d_all_signatures %>% 
    filter(value!=0) %>% 
    group_by(signature,tool,context) %>% 
    count() %>% 
    mutate(signature_freq_cohort=n/tot_samples) %>% 
    filter(signature_freq_cohort>=0.6) %>%
    ggplot()+
    geom_segment(aes(x = signature, y = 0, 
                    xend = signature, yend = signature_freq_cohort,color=signature)) +
    geom_point(size = 3, aes(x=signature,y=signature_freq_cohort,color=signature)) +
    xlab("Signature") +
    ylab("% in the cohort") +
    scale_fill_manual(values=COSMIC_color_palette())+
    scale_color_manual(values=COSMIC_color_palette())+
    facet_wrap(~tool,scales = "free")+
    my_ggplot_theme()+
    theme(legend.position = "none")+
    labs(title = "Signatures frequency in the cohort",caption = "Signatures present in more than 60% of the samples")


    final_report <- wrap_plots(list(all_cohort_exp,frq_cohort_plot,summary_tot_sig,
                                    exp_classes_all,
                                    cosmic_assign_sparse_sig_plot),
                               design = 'aaaaa\\naaaaa\\nbbbcc\\nddddd\\neeeee')




    # Save outputs
    write.table(d_all_signatures, file = paste0("$prefix", "_cohort_signature_exposures.txt"), sep = "\\t", row.names = FALSE, quote = FALSE)
    saveRDS(final_report, file = paste0("$prefix", "_cohort_signature_report.rds"))
    ggplot2::ggsave(plot = final_report, filename = paste0("$prefix", "_cohort_signature_report.pdf"), 
                    height = 12,width = 10)
    # Export versions
    f <- file("versions.yml", "w")
    ggplot2_version <- sessionInfo()\$otherPkgs\$ggplot2\$Version
    writeLines(paste0('"', "$task.process", '"', ":"), f)
    writeLines(paste("    ggplot2:", ggplot2_version), f)
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
