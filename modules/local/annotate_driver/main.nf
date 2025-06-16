process ANNOTATE_DRIVER {
    tag "$meta.id"
    label "process_single"
    label "error_retry"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://lvaleriani/cnaqc:version1.0' :
        'docker.io/lvaleriani/cnaqc:version1.0' }"

    input:
    tuple val(meta), path(rds), path(driver_list)

    output:
    tuple val(meta), path("*.rds"),     emit: rds
    path "versions.yml",                emit: versions

    script:
    def args    =   task.ext.args   ?: ''
    def prefix  =   task.ext.prefix ?: "${meta.id}"

    """
    #!/usr/bin/env Rscript
    library(dplyr)
    library(readr)
    library(tidyr)

    data = readRDS("$rds")
    SNV = data[["$meta.tumour_sample"]]
    SNV = SNV\$mutations

    drivers_table = readr::read_tsv(file = "$driver_list")

    tumour_type = "$meta.cancer_type"
    if(tumour_type %in% drivers_table\$TUMOUR_TYPE){
        drivers_table = drivers_table %>%
            dplyr::filter(TUMOUR_TYPE == tumour_type)
    } else {
        drivers_table = drivers_table %>%
            dplyr::mutate(TUMOUR_TYPE = "PANCANCER")
        tumour_type = 'PANCANCER'
    }

    drivers_table = drivers_table %>%
        dplyr::select(SYMBOL, TUMOUR_TYPE) %>%
        dplyr::distinct()

    x = SNV %>%
        dplyr::mutate(TUMOUR_TYPE = tumour_type) %>%
        dplyr::left_join(
            drivers_table,
            by = c('SYMBOL', 'TUMOUR_TYPE')
        ) %>%
        tidyr::separate(HGVSp, ':', into = c('s1', 's2'), remove=F) %>%
        dplyr::mutate(tmp_s2 = ifelse(is.na(s2), '', paste0('_', s2))) %>%
        dplyr::mutate(
            is_driver = (IMPACT %in% c('MODERATE', 'HIGH')),
            driver_label = paste0(SYMBOL, tmp_s2)
        ) %>%
        select(-tmp_s2) %>%
        mutate(is_driver = ifelse(is.na(is_driver), FALSE, is_driver))

    filter_x = x %>%
        distinct(chr, from, to, ref,  alt,  IMPACT, SYMBOL, Gene, is_driver, driver_label, .keep_all = T) %>%
        mutate(priority = ifelse(is_driver == TRUE, 1, 0)) %>%
        arrange(chr, from, to, desc(priority)) %>%
        distinct(chr, from, to, .keep_all = TRUE)

    data[["$meta.tumour_sample"]]\$mutations = filter_x
    saveRDS(object = data, file = paste0("$prefix", "_driver.rds"))

    # version export
    f <- file("versions.yml","w")
    readr_version <- sessionInfo()\$otherPkgs\$readr\$Version
    dplyr_version <- sessionInfo()\$otherPkgs\$dplyr\$Version
    tidyr_version <- sessionInfo()\$otherPkgs\$tidyr\$Version
    writeLines(paste0('"', "$task.process", '"', ":"), f)
    writeLines(paste("    readr:", readr_version), f)
    writeLines(paste("    dplyr:", dplyr_version), f)
    writeLines(paste("    tidyr:", tidyr_version), f)
    close(f)

    """
}
