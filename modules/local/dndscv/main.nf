//
// DNDSCV PROCESS
//
// dndscv module based on knitr prototype
// dndscv script (SINGLE SAMPLE ONLY)
// TODO multisamples
// 1. load input rds get mutations and sample id  
// 2. create dndscv input
// 3. Get list of drivers from file
// 4. Use custom reference

process DNDSCV {
  debug true
  tag "$meta.id"
  container "${workflow.containerEngine == 'singularity' ? 'docker://tucano/dndscv:latest' : 'tucano/dndscv:latest'}"


  input:
    
    tuple val(meta), path(snv_rds), path(driver_list), path(reference)
  
  output:

    tuple val(meta), path("*_dnds.rds"), emit: dnds_rds
  
  script:

    def args                                = task.ext.args    
    def prefix                              = task.ext.prefix                                       ?: "${meta.id}"  
    
    """
    #!/usr/bin/env Rscript

    library(stringr)
    library(dndscv)
    library(dplyr)
    library(readr)

    input <- readRDS("${snv_rds}")
    reference <- "${reference}"

    load(reference)
    
    sample_id <- input\$`${meta.tumour_sample}`\$sample
    mutations <- input\$`${meta.tumour_sample}`\$mutations
    
    # check if reference and mutations start with chr
    chr_ref <- as.character(RefCDS[[1]]\$chr)
    reference_with_chr <- startsWith(chr_ref, "chr")
    chr_mut <- as.character(mutations\$chr[1])
    mutations_with_chr <- startsWith(chr_mut,"chr")
    
    # always remove chr in dndscv_input
    dndscv_input <- mutations |>
      dplyr::select(chr,from,ref,alt) |>
      dplyr::rename(pos="from") |>
      dplyr::mutate(sample=sample_id, .before = chr) |>
      dplyr::mutate(chr=str_replace(chr,"chr",""))

    # add chr in dndscv_input when needed
    if (reference_with_chr) {
      dndscv_input |> dplyr::mutate(chr=paste("chr",chr,sep=""))
    }
    
    driver_genes <- read_lines("${driver_list}")

    dndscv_result <- dndscv::dndscv(
      mutations = dndscv_input,
      outmats = TRUE,
      max_muts_per_gene_per_sample = Inf,
      max_coding_muts_per_sample = Inf,
      outp = 3,
      use_indel_sites = TRUE,
      min_indels = 1,
      refdb=reference
    )

    sel <- left_join(dndscv_result\$sel_cv, dndscv_result\$sel_loc)
    annotation <- left_join(dndscv_result\$annotmuts, sel, by = c("gene"="gene_name"))
    annotation <- annotation |> 
      dplyr::select(!sampleID) |> 
      dplyr::rename(from="pos", alt="mut")

    # Add column known_driver based on driver list and potential_driver based on dndscv result
    annotation <- annotation |> 
      dplyr::mutate(
        known_driver=gene %in% driver_genes, 
        potential_driver = (qmis_cv <= 0.1 | qtrunc_cv <= 0.1 | qallsubs_cv <= 0.1)
      )
    

    # Here I must add chr when needed, as before 
    # I remove by default and reapply when needed
    mutations <- mutations |> dplyr::mutate(chr=str_replace(chr,"chr",""))
    annotation <- annotation |> dplyr::mutate(chr=str_replace(chr,"chr",""))
    if (mutations_with_chr) {
      mutations <- mutations |> dplyr::mutate(chr=paste("chr",chr,sep=""))
      annotation <- annotation |> dplyr::mutate(chr=paste("chr",chr,sep=""))
    }
    
    output_mutations <- left_join(mutations, annotation)

    output <- input
    output\$`${meta.tumour_sample}`\$mutations <- mutations
    
    saveRDS(output,"${prefix}_dnds.rds")
    """
}
