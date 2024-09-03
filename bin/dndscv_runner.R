#!/usr/bin/env Rscript

# Davide Rambaldi
# davide.rambaldi@gmail.com
# github.com/tucano

####################################################################################################################
# dndscv script (SINGLE SAMPLE ONLY)
# TODO multisamples
#
# Usage: ./bin/dndscv_runner.R -i ../tumourevo_datasets/VCF.rds -r ../tumourevo_datasets/RefCDS_human_GRCh37_dNdScv.0.1.0.rda -s Set6_42 -d ../tumourevo_datasets/driver_list.txt  -o test.rds
#
####################################################################################################################

library(optparse)
library(stringr)
library(readr)
library(dndscv)
library(dplyr)

option_list = list(
  make_option(
    c("-i", "--input"), 
    action="store", 
    default=NA, 
    type='character',
    help="input rds path"
  ),

  make_option(
    c("-o", "--output"), 
    action="store", 
    default=NA, 
    type='character',
    help="output rds path"
  ),

  make_option(
    c("-r", "--reference"), 
    action="store", 
    default=NA, 
    type='character',
    help="reference path"
  ),

  make_option(
    c("-s", "--sample"), 
    action="store", 
    default=NA, 
    type='character',
    help="sample id"
  ),

  make_option(
    c("-d", "--drivers"), 
    action="store", 
    default=NA, 
    type='character',
    help="drivers list"
  ),

  make_option(
    "--qmis_cv",
    default=0.1, 
    help="drivers list"
  ),
  
  make_option(
    "--qtrunc_cv",
    default=0.1, 
    help="drivers list"
  ),

  make_option(
    "--qallsubs_cv",
    default=0.1, 
    help="drivers list"
  )
)

opt = parse_args(OptionParser(option_list=option_list))

input <- readRDS(opt$i)
reference <- opt$r
load(reference)

sample_id <- opt$s
mutations <- input[[sample_id]]$mutations

# check if reference and mutations start with chr
chr_ref <- as.character(RefCDS[[1]][['chr']])
reference_with_chr <- startsWith(chr_ref, "chr")
chr_mut <- as.character(mutations[['chr']][1])
mutations_with_chr <- startsWith(chr_mut,"chr")

# always remove chr in dndscv_input
dndscv_input <- mutations |>
  dplyr::select(chr,from,ref,alt) |>
  dplyr::rename(pos="from") |>
  dplyr::mutate(sample=sample_id, .before = chr) |>
  dplyr::mutate(chr=str_replace(chr,"chr",""))

# add chr in dndscv_input when needed
if (reference_with_chr) {
  dndscv_input <- dndscv_input |> dplyr::mutate(chr=paste("chr",chr,sep=""))
}

driver_genes <- read_lines(opt$d)

dndscv_result <- dndscv::dndscv(
  mutations = dndscv_input,
  outmats = TRUE,
  max_muts_per_gene_per_sample = Inf,
  max_coding_muts_per_sample = Inf,
  outp = 3,
  use_indel_sites = TRUE,
  min_indels = 1,
  refdb=reference,
  cv=NULL
)

sel <- left_join(dndscv_result[['sel_cv']], dndscv_result[['sel_loc']])
annotation <- left_join(dndscv_result$annotmuts, sel, by = c("gene"="gene_name"))
annotation <- annotation |> 
  dplyr::select(!sampleID) |> 
  dplyr::rename(from="pos", alt="mut")

# Add column known_driver based on driver list and potential_driver based on dndscv result
annotation <- annotation |> 
  dplyr::mutate(
    known_driver=gene %in% driver_genes, 
    potential_driver = (qmis_cv <= opt$qmis_cv | qtrunc_cv <= opt$qtrunc_cv | qallsubs_cv <= opt$qallsubs_cv)
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
output[[sample_id]]$mutations <- mutations

# add globaldnds as a data.frame
output[[sample_id]]$globaldnds <- dndscv_result$globaldnds

# add theta 
# theta < 1 may reflect problems with the suitability of the dNdScv model for the dataset.
output[[sample_id]]$dndscv$theta <- dndscv_result$nbreg$theta
saveRDS(output,opt$o)