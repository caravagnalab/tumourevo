rm(list=ls())
.libPaths()
library(tidyverse)
library(CNAqc)
library(ComplexHeatmap)
library('ggalign')
library(patchwork)

# library(cowplot)
# library(gridExtra)
# library(grid)
# library(ggplotify)

# computing TMB per sample
data_path = '/orfeo/cephfs/scratch/cdslab/shared/tumourevo_crc_test/results_crc_cohort/driver_annotation/annotate_driver/CRCtest'

sample = list.files(data_path, full.names = T, recursive = T)[1]

sample = readRDS(sample)

sequenced_mb = 3.1e9/10^6

data = sample[[1]]$mutations %>%  # nb: need to find a way to access better the data --> but it is run on single samples so...
    separate(Consequence, into = "Consequence", sep = "&")

# compute TMB
tmb_stats = compute_tmb(data, seq_length = sequenced_mb) %>% 
  mutate(status = ifelse(TMB >= 10, 'Hyper-mutant (TMB >= 10)', 'TMB < 10'))

# plots

# plot number of mutations (or TMB) per chromosome
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

# add annotation of the drivers on the vaf

vaf_chr = data %>% 
  filter(chr %in% paste0('chr', c(seq(1:22), 'X', 'Y'))) %>% 
  mutate(chr = factor(chr, levels = (paste0('chr', c(seq(1:22), 'X', 'Y'))))) %>% 
  # mutate(driver_label = ifelse(is_driver == TRUE, driver_label, '')) %>% 
  ggplot(aes(
    VAF, 
    # label = driver_label
  )) + 
  geom_histogram(fill = '#B4D3D9') + 
  # ggrepel::geom_text_repel(
  #   # aes(label = driver_label, x = VAF),
  #   y = 10
  #   ) + 
  # ggrepel::geom_text_repel(y = .7) + 
  theme_bw() + 
  facet_wrap(~chr, scales = 'free_y', ncol = 6)

p_consequence = data %>% 
  separate(Consequence, into = "Consequence", sep = "&") %>% 
  ggplot(aes(fill = Consequence, 
             x = IMPACT)
             # y = n, 
             # fill = reorder(Consequence, -n))
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
  # guides(fill = guide_legend(ncol = 2))

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

# drivers in the sample --> w/ mut type?

# drivers heatmap
drivers = data %>% 
  filter(is_driver)

drivers_matrix = drivers %>% 
  mutate(VAF_class = ifelse(VAF <= .1, 'VAF < 0.1', 'VAF >= 0.1')) %>% 
  dplyr::select(driver_label, Consequence, VAF_class) %>% 
  mutate(mut = Consequence, 
         type = paste(Consequence, VAF_class, sep = ',')
         ) %>% 
  dplyr::select(-c(Consequence, VAF_class)) %>% 
  pivot_wider(values_from = type, #c(Consequence, VAF_class),
              names_from = driver_label) %>% 
  tibble::column_to_rownames('mut') %>% 
  as.matrix()

cols = setNames(
  nm = unique(drivers$Consequence), 
  object = RColorBrewer::brewer.pal('Spectral', n = length(unique(drivers$Consequence)))
)


pch_vaf = c(
  'VAF < 0.1' = 8, 
  'VAF >= 0.1' = ''
)

alter_fun = lapply(unique(drivers$Consequence), function(x) {
  alter_graphic("rect", width = 0.95, height = 0.95, fill = cols[x])
})
names(alter_fun) = unique(drivers$Consequence)

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

# ht_p = draw(ht, heatmap_legend_side = 'bottom')

ht_p <- draw(ht, heatmap_legend_side = "bottom")

design = '
AAAB
CCCD
'

page1 = wrap_plots(list(
  vaf_chr, p_chr, p_consequence, p_mut_type
), design = design)

pdf('test.pdf', width = 10, height = 10)
page1
ht_p
dev.off()




