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

# load libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(ComplexHeatmap)
library(patchwork)

# colors and utilities functions
compute_tmb = function(x, seq_length) {
  
  coding_muts = x %>% 
    separate(Consequence, into = "Consequence", sep = "&") %>% 
    filter(Consequence != 'synonymous_variant') %>% 
    filter(VAF > 0)
  
  coding_tmb = coding_muts %>% 
    group_by(sample) %>% 
    count() %>% 
    mutate(TMB = n/seq_length) 
  
  return(coding_tmb)
  
}

mut_cols = c('SNV' = '#7FBC41', 'Indel' = '#DE77AE')

# consequences
# colors of mutations
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
          'regulatory_region_variant', 'intergenic_variant', 'sequence_variant'), 
  object = c("wheat4", "darkseagreen2", "gold2", "mistyrose4", 'coral', "bisque", "mediumpurple", 
             "bisque3", "burlywood3", "blue3", "seashell4", "lightblue", "tan4", "orchid4", "darkorange", 
             "indianred", "seashell2", "plum3", "thistle2", "skyblue4", "red1", "darkolivegreen3", 
             "green", "tomato4", "turquoise4", "greenyellow", "cyan3", "slateblue3", "lightblue3", 
             "tomato", "sandybrown", "blue", "violetred4", "yellowgreen", "lightskyblue1", "blue2", 
             "salmon4", "darkseagreen1", "palegreen", "plum", "powderblue")
)

#computing TMB per sample
data = readRDS("$snv_rds") %>%
  purrr::pluck("$tumour_sample", "mutations") %>%
  dplyr::mutate(mutation_id = paste(chr,from,to,ref,alt,sep = ':'))

print(opt)

sequenced_mb = as.numeric(opt[['sequenced_mb']])/10^6

# compute TMB
t_type = data %>% 
  filter(!is.na(TUMOUR_TYPE)) %>%
  pull(TUMOUR_TYPE) %>% 
  unique
  
tmb_stats = compute_tmb(data, seq_length = sequenced_mb) %>% 
  mutate(status = ifelse(TMB >= 10, 'Hyper-mutant (TMB >= 10)', 'TMB < 10')) %>% 
  mutate(TUMOUR_TYPE = t_type)

# plots
# plot number of mutations per chromosome
p_chr = data %>% 
  filter(chr %in% paste0('chr', c(seq(1:22), 'X', 'Y'))) %>% 
  mutate(chr = factor(chr, levels = rev(paste0('chr', c(seq(1:22), 'X', 'Y'))))) %>% 
  ggplot(aes(
    chr
  )) + 
  geom_bar(stat = 'count', fill = '#B4D3D9') + 
  theme_bw() + 
  coord_flip() + 
  labs(x = 'Chromosome', 
       y = 'Number of mutations')

# plot vaf per chromosome
vaf_chr = data %>% 
  filter(chr %in% paste0('chr', c(seq(1:22), 'X', 'Y'))) %>% 
  mutate(chr = factor(chr, levels = (paste0('chr', c(seq(1:22), 'X', 'Y'))))) %>% 
  ggplot(aes(
    VAF
  )) + 
  geom_histogram(fill = '#B4D3D9') + 
  theme_bw() + 
  facet_wrap(~chr, scales = 'free_y', ncol = 6)

# plot type of mutation consequences
p_consequence = data %>% 
  separate(Consequence, into = "Consequence", sep = "&") %>% 
  ggplot(aes(fill = Consequence, 
             x = IMPACT)
         ) + 
  geom_bar(stat = 'count', position = 'stack') + 
  theme_bw() + 
  labs(x = '', 
       fill = 'Mutation effect', 
       y = 'Number of mutations', 
       fill = '')+
  theme(axis.text.y = element_blank(), 
        axis.ticks.length.y = unit(0, 'mm'), 
        legend.position = 'bottom'
        ) + 
  guides(fill = guide_legend(ncol = 6, title = '')) +
  facet_wrap(~IMPACT, scales = 'free', ncol = 1, strip.position = 'left') + # change the palette!
  coord_flip() + 
  scale_fill_manual(values = consequences_colors)

# plot the number of indels and snvs 
mut_type = data %>% 
  mutate(VARIANT_CLASS = case_when(
    VARIANT_CLASS %in% c('deletion', 'insertion') ~ 'Indel', 
    VARIANT_CLASS == 'substitution' ~ 'SNV', 
    .default = VARIANT_CLASS
    
  )) %>% 
  group_by(VARIANT_CLASS, sample) %>% 
  count()

# distribution of indels and snvs
p_mut_type = mut_type %>% 
  ggplot(aes(x = sample, 
             y = n, 
             fill = reorder(VARIANT_CLASS, -n))) + 
  geom_bar(stat = 'identity', position = 'dodge') + 
  theme_bw() + 
  labs(y = 'Number of mutations', 
       fill = '', 
       x = ' ') +
  scale_fill_manual(values = mut_cols) +
  theme(axis.text.x = element_blank(), 
        axis.ticks.length.x = unit(0, 'mm'), 
        legend.position = 'bottom'
        )

# drivers oncoprint (per sample)
drivers = data %>% 
  filter(is_driver)

drivers_matrix = drivers %>% 
  mutate(VAF_class = ifelse(VAF <= .1, 'VAF < 0.1', 'VAF >= 0.1')) %>% 
  dplyr::select(driver_label, Consequence, VAF_class) %>% 
  mutate(mut = Consequence, 
         type = paste(Consequence, VAF_class, sep = ',')
         ) %>% 
  dplyr::select(-c(Consequence, VAF_class)) %>% 
  pivot_wider(values_from = type, 
              names_from = driver_label) %>% 
  tibble::column_to_rownames('mut') %>% 
  as.matrix()

cols = setNames(
  nm = unique(drivers\$Consequence), 
  object = RColorBrewer::brewer.pal('Spectral', n = length(unique(drivers\$Consequence)))
)

pch_vaf = c(
  'VAF < 0.1' = 8, 
  'VAF >= 0.1' = ''
)

alter_fun = lapply(unique(drivers\$Consequence), function(x) {
  alter_graphic("rect", width = 0.95, height = 0.95, fill = cols[x])
})
names(alter_fun) = unique(drivers\$Consequence)

vaf_fun = lapply(pch_vaf, function(s) {
  function(x, y, w, h) 
    grid.points(x,y,w,h, pch = as.numeric(unname(s)), size = unit(3, "mm"))
})

alter_fun = c(alter_fun, 
              vaf_fun,
              list('background' = alter_graphic("rect", width = 0.95, height = 0.95, fill = 'gainsboro')))

ht = oncoPrint(drivers_matrix, 
          alter_fun = alter_fun, 
          col = cols, 
          show_row_names = F, 
          show_column_names = T, 
          top_annotation = NULL, 
          row_title = 'Consequence alterations', 
          column_title = 'Driver mutations', 
          heatmap_legend_param = list(ncol = 2), 
          name = 'Alteration type'
          )

# assembly everything together
design = '
AAAB
CCCD
'

page1 = wrap_plots(list(
  vaf_chr, p_chr, p_consequence, p_mut_type
), design = design)

ggplot2::ggsave(plot = page1, paste0(opt[['prefix']], '_mutations_report.pdf'), width = 260, height = 297, units="mm", dpi = 200)

# draw the oncoprint 

pdf(paste0(opt[['prefix']], '_driver_oncoprint.pdf'), width = 15, height = 8)
draw(ht, heatmap_legend_side = "bottom")
dev.off()

# save all the rds plots
saveRDS(object = tmb_stats, file = paste0(opt[['prefix']], '_tmb.rds'))
saveRDS(object = vaf_chr, file = paste0(opt[['prefix']], '_vaf_chr_plot.rds'))
saveRDS(object = p_chr, file = paste0(opt[['prefix']], '_chr_mut.rds'))
saveRDS(object = p_consequence, file = paste0(opt[['prefix']], '_consequence_mut_plot.rds'))
saveRDS(object = p_mut_type, file = paste0(opt[['prefix']], '_mut_type_plot.rds'))
saveRDS(object = ht, file = paste0(opt[['prefix']], '_driver_oncoprint.rds'))

# version export
f <- file("versions.yml","w")
dplyr_version <- sessionInfo()\$otherPkgs\$dplyr\$Version
ggplot2_version <- sessionInfo()\$otherPkgs\$ggplot2\$Version
tidyr_version <- sessionInfo()\$otherPkgs\$tidyr\$Version
complexheatmap_version <- sessionInfo()\$otherPkgs\$ComplexHeatmap\$Version
patchwork_version <- sessionInfo()\$otherPkgs\$patchwork\$Version
writeLines(paste0('"', "$task.process", '"', ":"), f)
writeLines(paste("    dplyr:", dplyr_version), f)
writeLines(paste("    tidyr:", tidyr_version), f)
writeLines(paste("    ggplot2:", ggplot2_version), f)
writeLines(paste("    ComplexHeatmap:", complexheatmap_version), f)
writeLines(paste("    patchwork:", patchwork_version), f)
close(f)