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
  prefix = ifelse('null' == 'null', 'MSeq', 'null'), 
  sequenced_mb = 3.1e9
)
opt_types <- lapply(opt, class)

# Apply parameter overrides

args_opt <- parse_args('--sequenced_mb 3.1e9')
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
library(maftools)
library(cowplot)
library(grid)


# colors and other setup stuff -----
cna_cols = CNAqc:::get_karyotypes_colors(c('0:0', '1:0', '1:1', '2:0', '2:1', '2:2'))
cna_cols = c(cna_cols, 'Other (Complex CNA)' = 'grey20')
names(cna_cols) = gsub(':', '-', names(cna_cols))

pch_mut = setNames(nm = c('Wild-type', 'Mutated'), 
                   c(NA, 8))

qc_colors = setNames(nm = c(TRUE, FALSE), 
  c('forestgreen', 'indianred3'))

tmb_cols = setNames(
  nm = c('Hyper-mutant (TMB >= 10)', 'TMB < 10'), 
  c('#BAABDA', '#D77FA1')
)

fga_cols = setNames(
  nm = c('FGA >= 0.2', 'FGA < 0.2'), 
  c('#FF7444', '#F2A65A')
)

purity_cols = setNames(
  nm = c('Purity <= 0.3', '0.3 < Purity <= 0.6', 'Purity > 0.6'), 
  object = c('#B4DEBD', '#91C4C3', '#80A1BA')
)

create_annotation = function(x, ann_colors, position) {
  
  if(!is.null(ann_colors)){
    ComplexHeatmap::HeatmapAnnotation(df = x, 
                                      col = ann_colors, 
                                      which = position, 
                                      show_annotation_name = F, 
                                      annotation_legend_param = list(nrow = 2, ncol = 3, width = 12, by_row = T))
  } else {
    ComplexHeatmap::HeatmapAnnotation(df = x, 
                                      which = position, 
                                      show_annotation_name = F, 
                                      annotation_legend_param = list(nrow = 2, ncol = 3, width = 12, by_row = T))
    }
}

compute_fga = function(cnaqc, sequenced) {
  fga = cnaqc\$basepairs_by_karyotype %>% 
    filter(karyotype != "1:1") %>%
    pull(n) %>% 
    sum()/sequenced
  return(fga)
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
tmb = lapply(strsplit('$tmb_rds', " ")[[1]], FUN = function(file){
  readRDS(file)
})
names(tmb) = tmb_patients
tmb = lapply(1:length(tmb), function(x) {
  tmb[[x]] %>% 
    mutate(Patient = names(tmb)[x])
}) %>% 
  bind_rows()

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

fga = lapply(original, function(or) {
  lapply(or, function(s) {
    tibble(
      sample = s\$sample,
      FGA = compute_fga(s, sequenced = opt[['sequenced_mb']]))
  }) %>% 
    bind_rows()
}) %>% 
  bind_rows()

ploidy_purity = lapply(original, function(or) {
  lapply(or, function(s) {
    tibble(
      sample = s\$sample,
      ploidy = s\$ploidy, 
      purity = s\$purity)
  }) %>% 
    bind_rows()
}) %>% 
  bind_rows()

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
  filter(VAF > 0) %>% 
  tidyr::separate(Consequence, sep = '&', into = 'Consequence', convert = T) %>% 
  dplyr::select(sample, SYMBOL, Consequence, karyotype, QC_PASS) %>% 
  group_by(sample, SYMBOL) %>% 
  mutate(Consequence = ifelse(n() > 1, 'multi-hit', Consequence)) %>% 
  distinct() %>% 
  mutate(karyotype = ifelse(karyotype %in% c('1:0', '1:1', '2:0', '2:1', '2:2'), karyotype, 'Other (Complex CNA)')) %>% 
  mutate(karyotype = gsub(':', '-', karyotype)) %>% 
  dplyr::select(-Consequence) %>% 
  distinct()

matrix_drivers = df %>% 
  # dplyr::filter(!is.na(Consequence)) %>%
  mutate(tb = paste(QC_PASS, karyotype, sep = ',')) %>%
  dplyr::select(-c(karyotype, QC_PASS)) %>% 
  pivot_wider(values_from = tb, names_from = SYMBOL) %>% 
  tibble::column_to_rownames('sample') %>% 
  t

# set alteration function for the mutation type
# dr_consequences = df\$Consequence %>% unique
# dr_consequences = dr_consequences[!is.na(dr_consequences)]
# cols = setNames(nm = dr_consequences, 
#                 object = consequences_colors[dr_consequences])
# 
# alter_fun = lapply(dr_consequences, function(x) {
#   alter_graphic("rect", width = 0.95, height = 0.95, fill = cols[x])
# })
# names(alter_fun) = dr_consequences
# 

# now add the karyotypes
karyos_data = df\$karyotype %>% unique 
cols_karyo = cna_cols[karyos_data]

cna_fun = lapply(names(cols_karyo), function(s) {
  alter_graphic("rect", width = 0.95, height = 0.95, fill = cols_karyo[s], alpha = .8)
})
names(cna_fun) = names(cols_karyo)

# adding background
cna_fun = c(cna_fun,
              list('background' = alter_graphic("rect", width = 0.95, height = 0.95, fill = 'gainsboro')))

qc_fun = lapply(names(qc_colors), function(qc) {
  alter_graphic("rect", width = 0.95, height = 0.95, col = qc_colors[qc], fill = NA, lty = 2, lwd = 3)
})
names(qc_fun) = names(qc_colors)

alter_fun = c(cna_fun, qc_fun)
# test_alter_fun(alter_fun)

# add annotations -- first one: patient

ann_top_data = tmb %>% 
  ungroup() %>% 
  full_join(., fga, by = 'sample') %>% 
  dplyr::select(Patient, sample, TMB, status, FGA) %>% 
  mutate(fga_status = ifelse(FGA >= .2, 'FGA >= 0.2', 'FGA < 0.2'))

samples_order = ann_top_data\$sample

ann_bottom_data = ploidy_purity[match(samples_order, ploidy_purity\$sample),] %>% 
  mutate(Purity = case_when(
    purity <= 0.3 ~ 'Purity <= 0.3', 
    (purity > 0.3 & purity <= 0.6) ~ '0.3 < Purity <= 0.6', 
    purity > 0.6 ~ 'Purity > 0.6'
  ))

bottom_ann = create_annotation(ann_bottom_data %>% dplyr::select(Purity),
                               ann_colors = list(Purity = purity_cols), 
                               position = 'column')
bottom_ann_ploidy <- HeatmapAnnotation(
  Ploidy = anno_barplot(
    ann_bottom_data\$ploidy,
    border = FALSE, 
    bar_width = 0.9
  ),
  annotation_name_side = "left"
)
 
bottom_ann = c(bottom_ann, bottom_ann_ploidy) 

# pt_cols = setNames(
#   nm = unique(ann_top_data\$Patient), 
#   object = RColorBrewer::brewer.pal(n = length(unique(ann_top_data\$Patient)), name = 'Pastel1')[1:length(unique(ann_top_data\$Patient))]
# )
# 
# ann_cols = list(
#  'Patient' = pt_cols
# )

top_ann = create_annotation(ann_top_data %>% dplyr::select(Patient), 
                            ann_colors = NULL,
                            position = 'column')

# tmb = list(
#   vec = setNames(nm = ann_top_data)
# )

top_ann_tmb <- HeatmapAnnotation(
  TMB = anno_barplot(
    setNames(object = ann_top_data\$TMB, nm = ann_top_data\$status),
    gp = gpar(fill = tmb_cols[ann_top_data\$status]),
    border = FALSE, 
    bar_width = 0.9
  ),
  annotation_name_side = "left"
)

top_ann_fga <- HeatmapAnnotation(
  FGA = anno_barplot(
    setNames(object = ann_top_data\$FGA, nm = ann_top_data\$fga_status),
    gp = gpar(fill = fga_cols[ann_top_data\$fga_status]),
    border = FALSE, 
    bar_width = 0.9
  ),
  annotation_name_side = "left"
)

top_ann = c(top_ann_fga, top_ann_tmb, top_ann)

lgd_tmb <- Legend(
  title = "TMB",
  at = names(tmb_cols),
  legend_gp = gpar(fill = tmb_cols)
)

lgs_fga = Legend(
  title = 'FGA', 
  at = names(fga_cols), 
  legend_gp = gpar(fill = fga_cols)
)

matrix_drivers = matrix_drivers[, samples_order]

ht = oncoPrint(matrix_drivers, 
               alter_fun = alter_fun, 
               col = cols_karyo, 
               show_row_names = T, 
               show_column_names = T, 
               top_annotation = top_ann, 
               bottom_annotation = bottom_ann, 
               #  row_title = 'Consequence alterations', 
               column_title = 'Driver mutations', 
               heatmap_legend_param = list(ncol = 2), 
               name = 'Copy number status and quality control'
)

pdf(paste0(opt[['prefix']],'_oncoprint.pdf'), width = 15, height = 12)
draw(ht, 
     heatmap_legend_side = 'bottom', 
     annotation_legend_side = 'bottom', 
     annotation_legend_list = c(list(lgd_tmb), list(lgs_fga)),
     merge_legend = TRUE)
dev.off()


# now create the TMB plot

# load tcga cohort data
tcga.cohort = system.file('extdata', 'tcga_cohort.txt.gz', package = 'maftools')
tcga.cohort = data.table::fread(file = tcga.cohort, sep = '\t', stringsAsFactors = FALSE)

tcga.cohort = tcga.cohort[,.(Tumor_Sample_Barcode, total, cohort)]
tcga.cohort$total = as.numeric(as.character(tcga.cohort$total))

# set each sample as a cohort name (or each patient?)
cohortName = unique(tmb$Patient)

tmb = tmb %>% 
  rename(cohort = Patient) %>%
  mutate(TCGA = 'Input')

samples_t_type = tmb$TUMOUR_TYPE %>% unique

tcga_t_types = as.data.frame(tcga.cohort) %>% 
  pull(cohort) %>% 
  unique

# associate tumour types of input samples with those from the tcga cohort and handle weird tumour types
tcga.cohort = lapply(samples_t_type, function(type) {
  
  if(!type %in% tcga_t_types) {
    
    tmb = tmb %>% 
      filter(TUMOUR_TYPE == type) %>% 
      mutate(TUMOUR_TYPE = 'PANCANCER')
    
    tcga.cohort = as.data.frame(tcga.cohort) %>% 
      rename(sample = Tumor_Sample_Barcode) %>% 
      rename(n = total) %>% 
      mutate(cohort = 'PANCANCER') %>% 
      mutate(TCGA = 'TCGA') %>% 
      # rename(Patient = cohort) %>% 
      mutate(TUMOUR_TYPE = cohort) %>% 
      bind_rows(., tmb) %>% 
      mutate(plot_total = n)
    
  } else {
    
    tmb = tmb %>% 
      filter(TUMOUR_TYPE == type)
    
    tcga.cohort = as.data.frame(tcga.cohort) %>% 
      rename(sample = Tumor_Sample_Barcode) %>% 
      rename(n = total) %>% 
      filter(cohort %in% samples_t_type) %>% 
      mutate(TCGA = 'TCGA') %>% 
      # rename(Patient = cohort) %>% 
      mutate(TUMOUR_TYPE = cohort) %>% 
      bind_rows(., tmb) %>% 
      mutate(plot_total = n) 
    
  }
    
}) %>% 
  bind_rows() %>% 
  distinct()

# #Median mutations
# tcga.cohort_median = tcga.cohort %>%
#   group_by(TUMOUR_TYPE) %>%
#   summarise(N = n(), Median_Mutations = median(plot_total)) %>%
#   arrange(if (decreasing) desc(Median_Mutations) else Median_Mutations) %>% 
#   rename(Cohort_size = N)

tcga.cohort = tcga.cohort %>% 
  mutate(TUMOUR_TYPE = factor(tcga.cohort$TUMOUR_TYPE, levels = tcga.cohort$TUMOUR_TYPE)) 

tcga.cohort = split(tcga.cohort, as.factor(tcga.cohort$cohort))
plot.dat = lapply(seq_len(length(tcga.cohort)), function(i){
  x = tcga.cohort[[i]]
  pos = rev(seq(i-1, i, length.out = nrow(x)))
  x %>% 
    arrange(desc(plot_total)) %>% 
    mutate(V1 = pos)  
}) %>% 
  bind_rows() %>% 
  mutate(cohort = factor(cohort)) %>% 
  mutate(facet_id = as.numeric(cohort)) 

# precompute ranges outside the mutate
tcga_global_range <- plot.dat %>% 
  filter(TCGA == "TCGA") %>% 
  summarise(min_v1 = min(V1), max_v1 = max(V1))

tcga_type_range <- plot.dat %>%
  filter(TCGA == "TCGA") %>%
  filter(TUMOUR_TYPE %in% samples_t_type) %>% 
  group_by(TUMOUR_TYPE) %>%
  summarise(min_v1 = min(V1), max_v1 = max(V1), .groups = "drop")

data = plot.dat %>%
  left_join(tcga_type_range, by = "TUMOUR_TYPE") %>%
  mutate(
    # if no TCGA match for this tumour type, fall back to global TCGA range
    min_v1 = if_else(is.na(min_v1), tcga_global_range$min_v1, min_v1),
    max_v1 = if_else(is.na(max_v1), tcga_global_range$max_v1, max_v1),
  ) %>%
  group_by(TUMOUR_TYPE) %>%
  mutate(
    V1_scaled = scales::rescale(V1, to = c(min_v1[1], max_v1[1])),
    V1_scaled = if_else(TCGA == "TCGA", V1, V1_scaled)
  ) %>%
  ungroup() %>%
  select(-min_v1, -max_v1) %>% 
  mutate(cc = ifelse(TCGA == 'TCGA', 'TCGA', as.character(cohort)))

bg <- data %>%
  dplyr::distinct(TUMOUR_TYPE, facet_id) %>%
  dplyr::mutate(col = ifelse(facet_id %% 2 == 0, "1", "2"))

med_df <- data %>%
  dplyr::group_by(TUMOUR_TYPE) %>%
  dplyr::summarise(med = median(log10(plot_total), na.rm = TRUE))

cex_opt = getOption('CNAqc_cex', default = 1)
plt_tmb <- data %>% 
  ggplot() +
  geom_rect(
    data = bg,
    aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = col),
    inherit.aes = FALSE,
    alpha = 0.1,
    show.legend = F
  ) +
  scale_fill_manual(values = bg_col) +
  geom_hline(
    data = med_df,
    aes(yintercept = med),
    color = "gray60",
    linewidth = 0.7
  ) + 
  geom_hline(
    yintercept = 0:6,
    linetype = "dashed",
    color = "grey",
    linewidth = 0.3
  ) + 
  geom_point(data = data  %>% filter(TCGA == 'TCGA'), aes(x = V1_scaled, y = log10(plot_total), col = cc), size =.4) +
  geom_point(data = data  %>% filter(TCGA != 'TCGA'), aes(x = V1_scaled,y = log10(plot_total), col = cc), size = 2) +
  # scale_color_manual('', values = col_point)+ 
  # scale_shape_manual('', values = c('Normal' = 16, 'WGD' =15, 'Hypermutant' = 17))+ 
  facet_grid(.~TUMOUR_TYPE, scales = 'free_x', switch = 'x') +
  ylab('TMB') +
  xlab('') +
  ggplot2::theme_light(base_size = 10 * cex_opt) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text.y = element_text(size = 8, margin = margin(l = 0, r = 0), colour = 'gray20'),
    panel.background = element_rect(fill = "transparent", colour = NA),     # removes grey background
    panel.grid.major = element_blank(),      # removes major grid lines
    panel.grid.minor = element_blank(),      # removes minor grid lines
    panel.spacing = unit(0, "mm"),                       # remove spacing between facets
    panel.border = element_blank(),
    plot.margin      = margin(0, 0, 0, 0), 
    plot.background = element_rect(color = "transparent", fill = NA, colour = NA),
    strip.background = element_blank(),
    strip.text.x = element_text(
      angle = 0,      # rotation
      vjust = 0.5,     # vertical alignment
      hjust = 0.5,      # horizontal alignment,
      size = 8, 
      margin = margin(t = 0.1, b = 0.1), 
      colour = 'gray20'
    )
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 2), title = ''),
    fill  = guide_legend(override.aes = list(size = 2))
    # shape = guide_legend(override.aes = list(size = 2))
  )  

ht = grid.grabExpr(draw(ht, 
          heatmap_legend_side = 'bottom', 
          annotation_legend_side = 'bottom', 
          annotation_legend_list = c(list(lgd_tmb), list(lgs_fga)),
          merge_legend = TRUE), 
          width  = 7,   # inches — pin this to avoid layout fighting
          height = 7)

# --- Combine the two plots ---

combined <- plot_grid(
  ht, plt_tmb,
  ncol        = 2,
  labels      = c("A", "B"),
  rel_widths  = c(1.6, 1)   # give oncoprint more room
)

ggsave(plot = plt_tmb, filename = paste0(opt[['prefix']],'_tmb.pdf'), bg = 'white', width = 6, height = 8, units = 'in')
ggsave(filename = paste0(opt[['prefix']],'_oncoprint_tmb.pdf'), bg = 'white', width = 17, height = 8, units = 'in')


# version export
f <- file("versions.yml","w")
dplyr_version <- sessionInfo()\$otherPkgs\$dplyr\$Version
ggplot2_version <- sessionInfo()\$otherPkgs\$ggplot2\$Version
tidyr_version <- sessionInfo()\$otherPkgs\$tidyr\$Version
complexheatmap_version <- sessionInfo()\$otherPkgs\$ComplexHeatmap\$Version
cnaqc_version <- sessionInfo()\$otherPkgs\$CNAqc\$Version
writeLines(paste0('"', "NFCORE_TUMOUREVO:TUMOUREVO:GENOME_INTERPRETER:COHORT_MUTATIONS", '"', ":"), f)
writeLines(paste("    dplyr:", dplyr_version), f)
writeLines(paste("    tidyr:", tidyr_version), f)
writeLines(paste("    ggplot2:", ggplot2_version), f)
writeLines(paste("    ComplexHeatmap:", complexheatmap_version), f)
writeLines(paste("    CNAqc:", cnaqc_version), f)
close(f)
