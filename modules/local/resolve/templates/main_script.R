#!/usr/bin/env Rscript

pkgs <- c("RESOLVE", "dplyr", "tidyr", "tibble", "purrr", "stringr", "ggplot2", "patchwork")
sapply(pkgs, require, character.only = TRUE)

parse_args <- function(x) {
  # Remove brackets
  x = gsub("\\\\[","",x)
  x = gsub("\\\\]","",x)
 
  # Split into key:value pairs
  args_list = unlist(strsplit(x, ", (?=[^)]*(?:\\\\(|\$))", perl=TRUE))

  # Convert to tibble and split key:value
  tibble(arg = args_list) %>%
    separate(arg, into = c("key", "value"), sep = ":", extra = "merge", fill = "right") %>%
    dplyr::mutate(across(everything(), ~str_trim(.))) %>%
    dplyr::filter(!is.na(key) & !is.na(value)) %>%
    tibble::deframe()
}


opt = list(
    prefix = ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix'),
    genome = "NULL",
    predefined_dbs_mbs = "FALSE",
    K = "1:10",
    nmf_runs = "100",
    num_processes = "all",
    cross_validation_entries = "0.01",
    cross_validation_repetitions = "50",
    cross_validation_iterations = "5",
    seed = "NULL"
)

args_opt = parse_args('$task.ext.args')
for ( ao in names(args_opt)) opt[[ao]] = args_opt[[ao]]


# Script #

num_procs_string <- opt[["num_processes"]]

if (is.null(num_procs_string)) {
    stop("Missing required option: num_processes")
} else if (num_procs_string == "all") {
    n_procs <- Inf
} else {
    n_procs <- as.integer(opt[["num_processes"]])
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

# Generate the patient vs mutation count matrix from mutation data
# Load reference human-genome specification.
# The user must select, among the available choices, the reference genome consistent with the mutation dataset.

load_genome = function(genome, input_data) {
    if (genome == "GRCh37") {
        library(BSgenome.Hsapiens.1000genomes.hs37d5)
        bsg = BSgenome.Hsapiens.1000genomes.hs37d5::hs37d5
        input_data <- input_data %>% mutate(chrom = str_remove(chrom,"^chr"))

    } else if (genome == "GRCh38") {
        library(BSgenome.Hsapiens.UCSC.hg38)
        bsg = BSgenome.Hsapiens.UCSC.hg38

        # Leave 'chrom' unchanged for GRCh38
    }
    return(list(bsg = bsg, input_data = input_data))
}

data_list = load_genome(opt[["genome"]], input_data)
bsg = data_list[["bsg"]]
input_data = data_list[["input_data"]]


split_multibase_to_snvs <- function(df) {
  df %>%
    # Process each row
    mutate(row_id = row_number()) %>%
    group_by(row_id) %>%
    group_map(~ {
      row <- .x

      ref_seq <- unlist(strsplit(row[["ref"]], split = ""))
      alt_seq <- unlist(strsplit(row[["alt"]], split = ""))

      # Sanity check: lengths must match for substitutions
      if (length(ref_seq) != length(alt_seq)) {
        return(row) # leave indels unchanged
      }

      # Find positions where bases differ
      diff_pos <- which(ref_seq != alt_seq)

      # Create one SNV row per differing position
      snv_rows <- map_dfr(diff_pos, function(pos) {
        tibble(
          sample = row[["sample"]],
          chrom  = row[["chrom"]],
          start  = row[["start"]] + (pos - 1),
          end    = row[["start"]] + (pos - 1),
          ref    = ref_seq[pos],
          alt    = alt_seq[pos]
        )
      })

      snv_rows
    }) %>%
    bind_rows()
}


snvs <- split_multibase_to_snvs(input_data)
mut_counts_mnv <- RESOLVE::getMNVCounts(data = snvs,
					 predefined_dbs_mbs = opt[["predefined_dbs_mbs"]])

mut_counts_sbs <- RESOLVE::getSBSCounts(data = input_data,
					reference = bsg)

indels <- input_data[nchar(input_data$ref) != nchar(input_data$alt), ]
mut_counts_id <- RESOLVE::getIDCounts(data = indels,
				      reference = bsg)


mut_counts_list <- list(
  SBS = mut_counts_sbs,
  MNV = mut_counts_mnv,
  ID = mut_counts_id
)

saveRDS(object = mut_counts_list, file = paste0(opt[["prefix"]], "_mut_counts.rds"))


### De-novo fit | Cross-validation ###

data(background)

run_signature_analysis <- function(count_matrix, background) {
  if (is.null(count_matrix) || nrow(count_matrix) == 0) {
    message("Skipping — no data in this matrix.")
    return(NULL)
  }

  # De-novo fit
  res_denovo <- RESOLVE::signaturesDecomposition(
						 x = count_matrix,
						 K = eval(parse(text=opt[["K"]])),
						 background_signature = background,
						 nmf_runs = as.integer(opt[["nmf_runs"]]),
						 num_processes = n_procs
                                                 )
  # Assignment
  res_assign <- signaturesAssignment(
				     x = count_matrix, 
				     beta = res_denovo[["beta"]][[1]])

  # Cross-validation
  res_cv <- RESOLVE::signaturesCV(
				  x = count_matrix,
				  beta = res_denovo[["beta"]],
				  cross_validation_entries = as.numeric(opt[["cross_validation_entries"]]),
				  cross_validation_iterations = as.integer(opt[["cross_validation_iterations"]]),
				  cross_validation_repetitions = as.integer(opt[["cross_validation_repetitions"]]),
				  num_processes = n_procs
                                 ) 

  list(denovo = res_denovo, assignment = res_assign, cv = res_cv)
}

# Apply to all available matrices
fit_results <- lapply(mut_counts_list, run_signature_analysis, background = background)

saveRDS(object = fit_results, file = paste0(opt[["prefix"]], "_fit_results.rds"))




