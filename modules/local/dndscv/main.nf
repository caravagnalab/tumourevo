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
  container='file:///fast/cdslab/ebusca00/singularity/cdslab.sif'

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
    
    # check reference chr mode
    reference_with_chr <- startsWith(RefCDS[[1]]\$chr, "chr")

    sample_id <- input\$`${meta.tumour_sample}`\$sample
    mutations <- input\$`${meta.tumour_sample}`\$mutations
    
    dndscv_input <- mutations |>
      dplyr::select(chr,from,ref,alt) |>
      dplyr::rename(pos="from") |>
      dplyr::mutate(sample=sample_id, .before = chr) |>
      dplyr::mutate(chr=str_replace(chr,"chr",""))

    # replace in dndscv_input and read if
    if (reference_with_chr) {
      dndscv_input |> dplyr::mutate(chr=paste("chr",chr,sep=""))
    }
    
    driver_genes <- read_lines("${driver_list}")

    # remove from driver genes the genes not in RefCDS
    refcds_genes = unlist(lapply(RefCDS, function(x){ x\$gene_name }), use.names=FALSE)
    driver_genes <- driver_genes[driver_genes %in% refcds_genes]

    dndscv_result <- dndscv::dndscv(
      mutations = dndscv_input,
      outmats = TRUE,
      max_muts_per_gene_per_sample = Inf,
      max_coding_muts_per_sample = Inf,
      outp = 3,
      use_indel_sites = TRUE,
      min_indels = 1,
      gene_list=driver_genes,
      refdb=reference
    )

    sel <- left_join(dndscv_result\$sel_cv, dndscv_result\$sel_loc)
    annotation <- left_join(dndscv_result\$annotmuts, sel, by = c("gene"="gene_name"))
    annotation <- annotation |> 
      dplyr::select(!sampleID) |> 
      dplyr::rename(from="pos", alt="mut")

    # TODO here I must add chr when needed
    output_mutations <- left_join(mutations, annotation)

    output <- input
    output\$`${meta.tumour_sample}`\$mutations <- mutations
    
    saveRDS(output,"${prefix}_dnds.rds")
    """
}
