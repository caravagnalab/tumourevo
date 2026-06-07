#!/usr/bin/env Rscript

library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(patchwork)

parse_args = function(x) {
  x = gsub("\\\\[","",x)
  x = gsub("\\\\]","",x)
  # giving errors when we have lists like c(xxx, xxx) since it will separate it
  # args_list = unlist(strsplit(x, ', ')[[1]])
  args_list = unlist(strsplit(x, ", (?=[^)]*(?:\\\\(|\$))", perl=TRUE))
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


get_cluster_colors <- function(cluster_names) {
  n <- length(cluster_names)

  color_pool <- unique(c(
    brewer.pal(8, "Dark2"),
    brewer.pal(12, "Paired")
  ))

  if (n > length(color_pool)) {
    color_pool <- colorRampPalette(color_pool)(n)
  }

  setNames(color_pool[1:n], nm = cluster_names)
}

if ("$rds_mobster" != ''){
  data = strsplit("$rds_mobster", ' ')[[1]]
  plt_mobster <- list()
  for (d in data){
    rds = readRDS(d)
    table = rds[["data"]] %>% select(VAF, SYMBOL, is_driver, cluster, sample_id)
    color = get_cluster_colors(unique(table[["cluster"]]))
    color[['Tail']] = 'gainsboro'
    plt = table %>%
      ggplot() +
      geom_histogram(aes(x = VAF, fill = cluster), binwidth = 0.01, position="identity", alpha = .7) +
      theme_minimal() +
      scale_fill_manual('Cluster', values = color) +
      xlim(0,1) +
      ggtitle(label = paste0('MOBSTER ', unique(table[["sample_id"]])), subtitle = paste0(nrow(table), ' mutations')) +
      ylab('')

    plt = plt +
      ggrepel::geom_label_repel(
        data = table %>% filter(is_driver == TRUE),
        aes(
          x = VAF,
          y = 0,
          label = SYMBOL,
          colour = cluster,
        ),
        show.legend = F,
        inherit.aes = FALSE,
        size = 3,
        min.segment.length = 0,
        box.padding = 1) +
      scale_color_manual('Cluster', values = color)
    plt_mobster[[unique(table[["sample_id"]])]] <- plt
  }
} else {
  plt_mobster = NULL
}

if ("$pyclone_fit" != ''){
  pyclone_table = read.table("$pyclone_table", header = T, sep = '\\t') %>%
    mutate(chr = sub("^chr", "", chr)) %>%
    mutate(mutation_id = paste(patient_id, chr, from, alt, sep = ':')) %>%
    select(mutation_id, VAF, SYMBOL, Indiv, is_driver) %>%
    dplyr::rename(sample_id = Indiv)
  pyclone_fit = read.table("$pyclone_fit", header = T, sep = '\\t') %>%
    select(mutation_id, cluster_id)


  plt_pyclone <- list()

  data = left_join(pyclone_table, pyclone_fit) %>%
    filter(!is.na(sample_id)) %>%
    mutate(cluster_id = paste0('C', cluster_id)) %>%
    distinct()

  samples = unique(data[["sample_id"]])
  if (length(samples)>1){
    pairs <- combn(samples, 2, simplify = FALSE)
    color = get_cluster_colors(unique(data[["cluster_id"]]))

    for (i in 1:length(pairs)){
      p = pairs[[i]]
      s1 = p[[1]]
      s2 = p[[2]]
      tmp = data %>%
        filter(sample_id %in% p) %>%
        pivot_wider(names_from = sample_id,
                    values_from = VAF)

      plt =  tmp %>%
        ggplot(aes(x =.data[[s1]], y = .data[[s2]], color=cluster_id))+
        geom_point( alpha=0.2, size = .5)+
        xlim(0,1)+
        ylim(0,1)+
        xlab(paste0('VAF ', s1)) +
        ylab(paste0('VAF ', s2)) +
        theme_minimal() +
        scale_color_manual('Cluster', values = color) +
        ggtitle(label = paste0('PyClone-VI'), subtitle = paste0(nrow(tmp), ' mutations')) +
        guides(color = guide_legend(override.aes = list(size = 2, alpha = 1)))


      plt = plt + ggrepel::geom_label_repel(
        data = tmp %>% filter(is_driver == TRUE),
        aes(
          x = .data[[s1]],
          y = .data[[s2]],
          label = SYMBOL,
          colour = cluster_id,
        ),
        show.legend = F,
        inherit.aes = FALSE,
        size = 3,
        min.segment.length = 0,
        box.padding = 1)

      plt_pyclone[[paste(s1, s2, sep = '-')]] <- plt

    }
  } else {
    sample = samples[[1]]
    color = get_cluster_colors(unique(data[["cluster_id"]]))

    plt = data %>%
      filter(VAF > 0) %>%
      ggplot() +
      geom_histogram(aes(x = VAF, fill = cluster_id), binwidth = 0.01, position="identity", alpha = .7) +
      theme_minimal() +
      scale_fill_manual('Cluster', values = color) +
      xlim(0,1) +
      ggtitle(label = paste0('PyClone-VI ', sample), subtitle = paste0(nrow(data), ' mutations')) +
      ylab('')

    plt = plt +
      ggrepel::geom_label_repel(
        data = data %>% filter(is_driver == TRUE) %>% distinct(),
        aes(
          x = VAF,
          y = 0,
          label = SYMBOL,
          colour = cluster_id,
        ),
        show.legend = F,
        inherit.aes = FALSE,
        size = 3,
        min.segment.length = 0,
        box.padding = 1) +
      scale_color_manual('Cluster', values = color)

    plt_pyclone[[sample]] <- plt
  }
} else {
  plt_pyclone = NULL
}

if ("$rds_viber" != ''){
  viber_fit = readRDS("$rds_viber")

  samples <- colnames(viber_fit[["data"]] %>%
                        select(starts_with("VAF"))) %>%
    str_remove(fixed("VAF."))

  data = viber_fit[["data"]] %>%
    select(starts_with("VAF"),gene, driver) %>%
    bind_cols(viber_fit[["labels"]])

  plt_viber <- list()

  if (length(samples)>1){
    pairs <- combn(samples, 2, simplify = FALSE)
    color = get_cluster_colors(unique(data[["cluster.Binomial"]]))

    for (i in 1:length(pairs)){
      p = pairs[[i]]
      s1 = p[[1]]
      s2 = p[[2]]
      vaf1 = paste0('VAF.', s1)
      vaf2 = paste0('VAF.', s2)

      plt =  data %>%
        ggplot(aes(x =.data[[vaf1]], y = .data[[vaf2]], color=cluster.Binomial))+
        geom_point( alpha=0.2, size = .5)+
        xlim(0,1)+
        ylim(0,1)+
        xlab(paste0('VAF ', s1)) +
        ylab(paste0('VAF ', s2)) +
        theme_minimal() +
        scale_color_manual('Cluster', values = color) +
        ggtitle(label = paste0('VIBER'), subtitle = paste0(nrow(data), ' mutations')) +
        guides(color = guide_legend(override.aes = list(size = 2, alpha = 1)))


      plt = plt + ggrepel::geom_label_repel(
        data = data %>% filter(driver == TRUE),
        aes(
          x = .data[[vaf1]],
          y = .data[[vaf2]],
          label = gene,
          colour = cluster.Binomial,
        ),
        show.legend = F,
        inherit.aes = FALSE,
        size = 3,
        min.segment.length = 0,
        box.padding = 1)

      plt_viber[[paste(s1, s2, sep = '-')]] <- plt

    }
  } else {
    sample = samples[[1]]
    vaf = paste0('VAF.', sample)
    color = get_cluster_colors(unique(data[["cluster.Binomial"]]))

    plt = data %>%
      filter(.data[[vaf]] > 0.01) %>%
      ggplot() +
      geom_histogram(aes(x = .data[[vaf]], fill = cluster.Binomial), binwidth = 0.01, position="identity", alpha = .7) +
      theme_minimal() +
      scale_fill_manual('Cluster', values = color) +
      xlim(0,1) +
      ggtitle(label = paste0('VIBER ', sample), subtitle = paste0(nrow(data), ' mutations')) +
      ylab('') +
      xlab('VAF')

    plt = plt +
      ggrepel::geom_label_repel(
        data = data %>% filter(driver == TRUE) %>% distinct(),
        aes(
          x = .data[[vaf]],
          y = 0,
          label = gene,
          colour = .data[["cluster.Binomial"]],
        ),
        show.legend = F,
        inherit.aes = FALSE,
        size = 3,
        min.segment.length = 0,
        box.padding = 1) +
      scale_color_manual('Cluster', values = color)

    plt_viber[[sample]] <- plt
  }
} else {
  plt_viber = NULL
}


if (length(plt_mobster) == 1){
  ggsave(plt_mobster, filename=paste0(opt[["prefix"]], "_deconvolution_mobster.pdf"), width = 4, height = 3, units = 'in')

} else if (length(plt_mobster) > 1) {
  hg = ceiling(length(plt_mobster)/2)
  mobster <- wrap_plots(plt_mobster, ncol = 2, nrow = hg)
  ggsave(mobster, filename=paste0(opt[["prefix"]], "_deconvolution_mobster.pdf"), width = 8, height = hg*3, units = 'in')

}

if (length(plt_pyclone) == 1){
  ggsave(plt_pyclone, filename=paste0(opt[["prefix"]], "_deconvolution_pyclonevi.png"), width = 4, height = 3, units = 'in')

} else if (length(plt_pyclone) > 1) {
  hg = ceiling(length(plt_pyclone)/2)
  pyclone <- wrap_plots(plt_pyclone, ncol = 2, nrow = hg)
  ggsave(pyclone, filename=paste0(opt[["prefix"]], "_deconvolution_pyclonevi.png"), width = 8, height = hg*3, units = 'in')

}

if (length(plt_viber) == 1){
  ggsave(plot = plt_viber, filename=paste0(opt[["prefix"]], "_deconvolution_viber.png"), width = 4, height = 3, units = 'in')

} else if (length(plt_viber) > 1){
  hg = ceiling(length(plt_viber)/2)
  viber <- wrap_plots(plt_viber, ncol = 2, nrow = hg)
  ggsave(plot = viber, filename=paste0(opt[["prefix"]], "_deconvolution_viber.png"), width = 8, height = hg*3, units = 'in')

}
