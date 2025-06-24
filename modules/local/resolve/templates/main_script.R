#!/usr/bin/env Rscript

pkgs <- c("RESOLVE", "dplyr", "tidyr", "tibble", "purrr", "stringr", "ggplot2", "patchwork")
sapply(pkgs, require, character.only = TRUE)


parse_args = function(x) {
    x = gsub("\\\\[","",x)
    x = gsub("\\\\]","",x)
    # giving errors when we have lists like c(xxx, xxx) since it will separate it
    # args_list = unlist(strsplit(x, ', ')[[1]])
    args_list = unlist(strsplit(x, ", (?=[^)]*(?:\\\\(|\$))", perl=TRUE))
    # args_vals = lapply(args_list, function(x) strsplit(x, split=":")[[1]])
    args_vals = lapply(args_list, function(x) {
        x_splt = strsplit(x, split=":")[[1]]
        c(x_splt[1],  paste(x_splt[2:length(x_splt)], collapse=":"))
    })

    # Ensure the option vectors are length 2 (key/ value) to catch empty ones
    args_vals = lapply(args_vals, function(z){ length(z) = 2; z})

    parsed_args = structure(lapply(args_vals, function(x) x[2]), names = lapply(args_vals, function(x) x[1]))
    parsed_args[! is.na(parsed_args)]
}

opt = list(
    prefix = ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix')
)
args_opt = parse_args('$task.ext.args')
for ( ao in names(args_opt)) opt[[ao]] = args_opt[[ao]]


# Script

n_procs = parse(text=opt[["num_processes"]])
if (n_procs == "all"){
    n_procs = as.double("Inf")
} else {
    n_procs = eval(n_procs)
}


patients_tsv = strsplit("$tsv_join", " ")[[1]]
tables = lapply(patients_tsv, FUN = function(p_table){
    read.delim(p_table, sep = "\\t", header=T) %>%
        mutate(across(everything(), as.character))
}
)
multisample_table = dplyr::bind_rows(tables)


#Extract input data information
input_data = multisample_table[,c("Indiv","chr","from","to","ref","alt")]
input_data = setNames(input_data, c("sample","chrom","start","end","ref","alt"))
input_data[["end"]] = input_data[["start"]]
input_data = input_data %>% mutate(start = as.integer(start), end = as.integer(end))

#Generate the patient vs mutation count matrix from mutation data
#Load reference human-genome specification.
#The user must select, among the available choices, the reference genome consistent with the mutation dataset.

load_genome = function(genome, input_data) {
    if (genome == "GRCh37") {
        library(BSgenome.Hsapiens.1000genomes.hs37d5)
        bsg = BSgenome.Hsapiens.1000genomes.hs37d5::hs37d5
        input_data[["chrom"]] = substring(input_data[["chrom"]], 4, 5)

    } else if (genome == "GRCh38") {
        library(BSgenome.Hsapiens.UCSC.hg38)
        bsg = BSgenome.Hsapiens.UCSC.hg38

        # Leave 'chrom' unchanged for GRCh38
    }
    return(list(bsg = bsg, input_data = input_data))
}

data_list = load_genome("$params.genome", input_data)
bsg = data_list[["bsg"]]
input_data = data_list[["input_data"]]
