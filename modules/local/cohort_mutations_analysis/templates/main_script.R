#!/usr/bin/env Rscript

# parse arguments
parse_args <- function(x){
  args_list <- unlist(strsplit(x, ' ?--')[[1]])[-1]
  args_vals <- lapply(args_list, function(x) scan(text=x, what='character', quiet = TRUE))

  # Ensure the option vectors are length 2 (key/ value) to catch empty ones
  args_vals <- lapply(args_vals, function(z){ length(z) <- 2; z})

  parsed_args <- structure(lapply(args_vals, function(x) x[2]), names = lapply(args_vals, function(x) x[1]))
  parsed_args[! is.na(parsed_args)]
}

# Set defaults and classes

opt <- list(
  prefix = ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix'), 
  sequenced_mb = 3.1e9
)
opt_types <- lapply(opt, class)

# Apply parameter overrides

args_opt <- parse_args('$task.ext.args')
for ( ao in names(args_opt)){
  if (! ao %in% names(opt)){
    stop(paste("Invalid option:", ao))
  }else{

    # Preserve classes from defaults where possible
    if (! is.null(opt[[ao]])){
      args_opt[[ao]] <- as(args_opt[[ao]], opt_types[[ao]])
    }
    opt[[ao]] <- args_opt[[ao]]
  }
}


library(ComplexHeatmap)
library(dplyr)
library(CNAqc)
library(tidyr)
library(ggplot2)

# colors for mutation types 

consequences_colors = setNames(
  nm = c('transcript_ablation', 'splice_acceptor_variant', 'splice_donor_variant', 'stop_gained',    
         'frameshift_variant', 'stop_lost', 'start_lost', 'transcript_amplification', 'feature_elongation',
         'feature_truncation', 'inframe_insertion', 'inframe_deletion', 'missense_variant', 
         'protein_altering_variant', 'splice_donor_5th_base_variant', 'splice_region_variant', 
         'splice_donor_region_variant', 'splice_polypyrimidine_tract_variant', 
         'incomplete_terminal_codon_variant', 'start_retained_variant', 'stop_retained_variant',     
         'synonymous_variant', 'coding_sequence_variant', 'mature_miRNA_variant', '5_prime_UTR_variant',
         '3_prime_UTR_variant', 'non_coding_transcript_exon_variant', 'intron_variant', 
         'NMD_transcript_variant', 'non_coding_transcript_variant', 'coding_transcript_variant', 
         'upstream_gene_variant', 'downstream_gene_variant', 'TFBS_ablation', 'TFBS_amplification',    
         'TF_binding_site_variant', 'regulatory_region_ablation', 'regulatory_region_amplification',     
         'regulatory_region_variant', 'intergenic_variant', 'sequence_variant', 'multi-hit'), 
  object = c("wheat4", "darkseagreen2", "gold2", "mistyrose4", 'coral', "bisque", "mediumpurple", 
             "bisque3", "burlywood3", "blue3", "seashell4", "lightblue", "tan4", "orchid4", "darkorange", 
             "indianred", "seashell2", "plum3", "thistle2", "skyblue4", "red1", "darkolivegreen3", 
             "green", "tomato4", "turquoise4", "greenyellow", "cyan3", "slateblue3", "lightblue3", 
             "tomato", "sandybrown", "blue", "violetred4", "yellowgreen", "lightskyblue1", "blue2", 
             "salmon4", "darkseagreen1", "palegreen", "plum", "powderblue", 'seagreen')
)

pch_cna = setNames(
  nm = c('1-0', '1-1', '2-1', '2-0', 'Other'),
  object = c(16, 8, 15, 17, 4)
)

create_annotation = function(x, ann_colors, position) {
  
  ComplexHeatmap::HeatmapAnnotation(df = x, 
                                    col = ann_colors, 
                                    which = position, 
                                    show_annotation_name = F, 
                                    annotation_legend_param = list(nrow = 2, ncol = 3, width = 12, by_row = T))
  
}


# loading patients names

cnaqc_patients = substr("$cnaqc_patients", 2, nchar("$cnaqc_patients")-1)
cnaqc_patients = strsplit(cnaqc_patients, ", ")[[1]]

print(cnaqc_patients)

tmb_patients = substr("$tmb_patients", 2, nchar("$tmb_patients")-1)
tmb_patients = strsplit(tmb_patients, ", ")[[1]]

print(tmb_patients)

# load data

cnaqc_list = lapply(strsplit("$join_cnaqc", " ")[[1]], FUN = function(file){
            readRDS(file)
        })
names(cnaqc_list) = cnaqc_patients

# get the tmb results 
tmb = lapply(strsplit("$tmb_rds", " ")[[1]], FUN = function(file){
            readRDS(file)
        })
names(tmb) = tmb_patients

original = lapply(cnaqc_list, function(x) {
  df = CNAqc::get_sample(x, CNAqc::get_sample_name(x), which_obj = 'original')
  
  tt = lapply(CNAqc::get_sample_name(x), function(sp) {
    df[[sp]]\$sample = sp
    return(df[[sp]])
  })
  
  names(tt) = CNAqc::get_sample_name(x)
  return(tt)

  # return(df)
})
print(original)

muts = lapply(original %>% names, function(pt) {
  driver_muts = lapply(original[[pt]], function(sp) {
    
    drivers = CNAqc::Mutations(sp) %>% 
      mutate(sample = CNAqc::get_sample_name(sp)) %>% 
      filter(is_driver)
    
    driver_ccf = CNAqc::CCF(sp) %>% 
      filter(is_driver) %>% 
      mutate(sample = CNAqc::get_sample_name(sp)) 
    
    full_join(drivers, driver_ccf)  
    
  }) %>% 
    bind_rows() %>% 
    mutate(patient = pt)
  
}) %>% 
  bind_rows()

drivers = muts\$SYMBOL %>% unique

drivers_cna_status = lapply(original %>% names, function(pt) {
  lapply(original[[pt]], function(sp) {
    CNAqc::CNA_gene(sp, genes = drivers) %>% 
      mutate(sample = CNAqc::get_sample_name(sp)) 
  }) %>% 
    bind_rows() %>% 
    mutate(patient = pt)
  # dplyr::mutate(sample = get_sample_name(or))
}) %>% 
  bind_rows() %>% 
  rename(SYMBOL = gene) %>% 
  dplyr::select(SYMBOL, karyotype, sample, patient)

df = full_join(muts, drivers_cna_status) %>% 
  tidyr::separate(Consequence, sep = '&', into = 'Consequence', convert = T) %>% 
  dplyr::select(sample, SYMBOL, Consequence, karyotype) %>% 
  group_by(sample, SYMBOL) %>% 
  mutate(Consequence = ifelse(n() > 1, 'multi-hit', Consequence)) %>% 
  distinct() %>% 
  mutate(karyotype = ifelse(karyotype %in% c('1:0', '1:1', '2:0', '2:1'), karyotype, 'Other')) %>% 
  mutate(karyotype = gsub(':', '-', karyotype))

matrix_drivers = df %>% 
  dplyr::filter(!is.na(Consequence)) %>%
  mutate(mut_effect = paste(Consequence, karyotype, sep = ',')) %>%
  dplyr::select(-c(karyotype, Consequence)) %>% 
  pivot_wider(values_from = mut_effect, names_from = SYMBOL) %>% 
  tibble::column_to_rownames('sample') %>% 
  t

# set alteration function for the mutation type
dr_consequences = df\$Consequence %>% unique
dr_consequences = dr_consequences[!is.na(dr_consequences)]
cols = setNames(nm = dr_consequences, 
                object = consequences_colors[dr_consequences])

alter_fun = lapply(dr_consequences, function(x) {
  alter_graphic("rect", width = 0.95, height = 0.95, fill = cols[x])
})
names(alter_fun) = dr_consequences

# adding background
alter_fun = c(alter_fun, 
              list('background' = alter_graphic("rect", width = 0.95, height = 0.95, fill = 'gainsboro')))

# now add the karyotypes
karyos_data = df\$karyotype %>% unique 
pch_cna = pch_cna[karyos_data]

cna_fun = lapply(pch_cna, function(s) {
  function(x, y, w, h)
    grid.points(x,y,w,h, pch = as.numeric(unname(s)), size = unit(3, "mm"))
})
cna_fun = c(cna_fun, "NA" = function(x, y, w, h) {
  grid.points(x,y,w,h, pch = "", size = unit(3, "mm"), gpar(col = NA, alpha = 0, fill = NA))
})

alter_fun = c(alter_fun, cna_fun)

# add annotations -- first one: patient
patients = muts %>% 
  dplyr::select(sample, patient) %>% 
  distinct() %>% 
  rename(Patient = patient) %>% 
  dplyr::select(Patient)

pt_cols = setNames(
  nm = unique(patients\$Patient), 
  object = RColorBrewer::brewer.pal(n = length(unique(patients\$Patient)), name = 'Pastel1')[1:length(unique(patients\$Patient))]
)

ann_cols = list(
  'Patient' = pt_cols
)

top_ann = create_annotation(patients, 
                  ann_colors = ann_cols,
                  position = 'column')

ht = oncoPrint(matrix_drivers, 
               alter_fun = alter_fun, 
               col = cols, 
               show_row_names = T, 
               show_column_names = T, 
               top_annotation = top_ann, 
              #  row_title = 'Consequence alterations', 
               column_title = 'Driver mutations', 
               heatmap_legend_param = list(ncol = 2), 
               name = 'Alteration type (SNV and CNA)'
)

pdf(paste0(opt[['prefix']],'_oncoprint.pdf'), width = 15, height = 12)
draw(ht, heatmap_legend_side = 'bottom', annotation_legend_side = 'bottom', merge_legend = TRUE)
dev.off()

# version export
f <- file("versions.yml","w")
dplyr_version <- sessionInfo()\$otherPkgs\$dplyr\$Version
ggplot2_version <- sessionInfo()\$otherPkgs\$ggplot2\$Version
tidyr_version <- sessionInfo()\$otherPkgs\$tidyr\$Version
complexheatmap_version <- sessionInfo()\$otherPkgs\$ComplexHeatmap\$Version
cnaqc_version <- sessionInfo()\$otherPkgs\$CNAqc\$Version
writeLines(paste0('"', "$task.process", '"', ":"), f)
writeLines(paste("    dplyr:", dplyr_version), f)
writeLines(paste("    tidyr:", tidyr_version), f)
writeLines(paste("    ggplot2:", ggplot2_version), f)
writeLines(paste("    ComplexHeatmap:", complexheatmap_version), f)
writeLines(paste("    CNAqc:", cnaqc_version), f)
close(f)