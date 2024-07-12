process VCF_PROCESSING {
    tag "$meta.id"

    container "docker://lvaleriani/cnaqc:dev1"

    input:
        tuple val(meta), path(vcf), path(tbi)
    
    output:
        tuple val(meta), path("*.rds"), emit: rds

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def prefix = task.ext.prefix ?: "${meta.id}"
        def filter_mutations  = args!='' && args.filter_mutations     ?  "$args.filter_mutations" : ""

    """
    #!/usr/bin/env Rscript 
    library(tidyverse)
    library(vcfR)
    
    source(paste0("$moduleDir", '/parser_vcf.R'))

    # Read vcf file
    vcf = vcfR::read.vcfR("$vcf")

    # Check from which caller the .vcf has been produced
    source = vcfR::queryMETA(vcf, element = 'source')[[1]]

    if (grepl(pattern = 'Mutect', x = source)){
        calls = parse_Mutect(vcf, sample_id = "$meta.tumour_sample", filter_mutations = as.logical("$filter_mutations"))
        
    } else if (grepl(pattern = 'strelka', x = source)){
        calls = parse_Strelka(vcf, sample_id = "$meta.tumour_sample", filter_mutations = as.logical("$filter_mutations"))
    
    } else if (grepl(pattern = 'Platypus', x = source)){
        calls = parse_Platypus(vcf, sample_id = "$meta.tumour_sample", filter_mutations = as.logical("$filter_mutations"))

    } else if (grepl(pattern = 'freeBayes', x = source)){
        calls = parse_Freebayes(vcf, sample_id = "$meta.tumour_sample", filter_mutations = as.logical("$filter_mutations"))

    } else {
        stop('Variant Caller not supported.')
    }

    saveRDS(object = calls, file = paste0("$prefix", ".rds"))
    """
}
