#!/usr/bin/env Rscript

parse_args = function(x) {
  x = gsub("\\\\[","",x) # nolint: indentation_linter.
  x = gsub("\\\\]","",x)
  # errors when we have lists like c(xxx, xxx) since it will separate it
  # args_list = unlist(strsplit(x, ', ')[[1]])
  args_list = unlist(strsplit(x, ", (?=[^)]*(?:\\\\(|\$))", perl=TRUE))
  args_vals = lapply(args_list, function(x) {
    x_splt = strsplit(x, split=":")[[1]] # nolint: indentation_linter.
    c(x_splt[1],  paste(x_splt[2:length(x_splt)], collapse=":"))
  })

  # Ensure the option vectors are length 2 (key/ value) to catch empty ones
  args_vals = lapply(args_vals, function(z){ length(z) = 2; z})

  parsed_args = structure(lapply(args_vals, function(x) x[2]),
                          names = lapply(args_vals, function(x) x[1]))
  parsed_args[! is.na(parsed_args)]
}

opt = list(
  prefix = ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix')
)
args_opt = parse_args('$task.ext.args')
for ( ao in names(args_opt)) opt[[ao]] = args_opt[[ao]]

# Script #####

mutation_tables = strsplit(gsub("\\[|\\]", "", "$mutation_tables"), ", ")[[1]]
results_sigprofiler = strsplit(gsub("\\[|\\]", "", "$results_sigprofiler"), ", ")[[1]]

library(tidyverse)
library(ggplot2)

signature_colors = c("#f1696bff", "#8fbd8cff", "#87c7d6ff", "#bac3deff",
                     "#d7bfd9ff", "#a8a2a1ff", "#cfadb3ff", "#3c609aff",
                     "#9a4564ff", "#fbcb5bff", "#c2b280ff", "#d47e2dff",
                     "#5f8676ff", "forestgreen", "orange", "brown4")

colors_cluster = c("indianred", "steelblue", "forestgreen", "goldenrod",
                   "darkorange3", "palevioletred", "mediumpurple", "cornsilk4",
                   "olivedrab3", "steelblue4", "indianred4", "aquamarine3",
                   "saddlebrown", "deeppink2", "cornflowerblue", "black") %>%
  setNames(paste0("C",0:15))


cosine_similarity = function(vec1, vec2) {
  sum(vec1 * vec2) / (sqrt(sum(vec1^2)) * sqrt(sum(vec2^2)))
}

compare_signatures = function(df1, df2) {
  df1_f = df1 %>% filter(Exposure>.1)
  df2_f = df2 %>% filter(Exposure>.1)

  # Check if signatures are identical
  sigs_match = setequal(df1_f$Signature, df2_f$Signature)
  ndiff = length(setdiff(df1_f$Signature, df2_f$Signature))
  ntot = length(c(df1_f$Signature, df2_f$Signature) %>% unique())

  comparison = full_join(df1 %>% select(Signature, Exposure, n_sig),
                         df2 %>% select(Signature, Exposure, n_sig),
                         by="Signature", suffix=c("_tmp", "_clonal")) %>%
    mutate(diff_exposure=Exposure_tmp-Exposure_clonal,
           diff_n_sig=n_sig_tmp-n_sig_clonal) %>%
    mutate(across(.cols=where(is.numeric),
                  .fns=function(i) replace_na(i, 0.)))

  comparison_data = full_join(df1, df2, by="Signature",
                              suffix=c("_tmp", "_clonal")) %>%
    select(Signature, Exposure_tmp, Exposure_clonal, n_sig_tmp, n_sig_clonal) %>%
    mutate(across(where(is.numeric), ~replace_na(., 0)))

  # Calculate similarities
  cos_sim_exposure = cosine_similarity(comparison_data$Exposure_tmp,
                                       comparison_data$Exposure_clonal)

  return(list("df"=comparison,
              "match"=sigs_match,
              "cs_exp"=cos_sim_exposure,
              "n_diff"=ndiff,
              "n_tot"=ntot))
}

# Signature interpretation ######

# patient_ids = list.dirs(mutation_tables, full.names=F, recursive=F)
# Get the patient ids from the mutation tables paths
patient_ids = lapply(mutation_tables, function(x) {
  tmp = strsplit(x, "/")[[1]]
  tmp[length(tmp)-1]
}) %>% unlist() %>% unique()

score_table = lapply(patient_ids, function(patient_id) {
  # mobster_files = list.files(file.path(mutation_tables, patient_id), recursive=T,
  #                            pattern="mobster", full.names=T)
  mobster_files = grep(mutation_tables, pattern="mobster", value=T)
  mutations_mobster = lapply(mobster_files, readr::read_tsv) %>%
    bind_rows() %>%
    mutate(patient_id=patient_id, tool="mobster",
           mutation_id=paste(chrom, pos_start, pos_end, ref, alt, sep=":")) %>%
    select(Project, patient_id, mutation_id, everything(), -chrom, -pos_start,
           -pos_end, -ref, -alt, -Type, -ID, -Genome, -mut_type) %>%
    rename(cluster_mobster=Sample)

  lapply(c("viber", "pyclonevi"), function(tool) {
    lapply(c("SBS", "ID"), function(sign_type) {
      # mutations_file = list.files(file.path(mutation_tables, patient_id), recursive=T, pattern=tool, full.names=TRUE)
      mutations_file = grep(mutation_tables, pattern=tool, value=T)
      if (length(mutations_file) == 0) return(tibble())
      mutations_tool = readr::read_tsv(mutations_file) %>%
        mutate(patient_id=patient_id, tool=tool,
               mutation_id=paste(chrom, pos_start, pos_end, ref, alt, sep=":")) %>%
        select(Project, patient_id, mutation_id, everything(), -chrom, -pos_start,
               -pos_end, -ref, -alt, -Type, -ID, -Genome, -mut_type) %>%
        rename(cluster_tool=Sample)

      never_tail_muts = mutations_mobster %>%
        group_by(mutation_id) %>% 
        summarise(never_tail=all(cluster_mobster != "Tail"))

      final_table_subclonal = mutations_tool %>%
        left_join(never_tail_muts) %>%
        group_by(cluster_tool) %>%
        reframe(n_never_tail=sum(never_tail, na.rm=T) / n(),
                is_driver=any(is_driver), is_clonal=all(is_clonal))

      # one dir per sample
      dir_path = list.dirs(grep(results_sigprofiler, pattern=tool, value=T), recursive=T) %>%
        keep(function(x) grepl(pattern=sign_type, x=x)) %>%
        keep(function(x) grepl(pattern="/Activities", x=x))

      if (length(dir_path) == 0) return(tibble())

      signatures = read.table(file=file.path(dir_path, "Assignment_Solution_Activities.txt"), header=T) %>%
        rename(cluster_tool=Samples) %>%
        pivot_longer(cols=starts_with(sign_type), names_to="Signature", values_to="Nmuts") %>%
        group_by(cluster_tool) %>%
        mutate(Ntot=sum(Nmuts), Exposure=Nmuts/Ntot, patient_id=patient_id) %>%
        ungroup() %>%
        left_join(mutations_tool %>% select(cluster_tool, is_clonal) %>% unique())

      clonal_sigs = signatures %>% filter(is_clonal, Nmuts>0) %>% rename(n_sig=Nmuts)
      other_sigs = signatures %>% filter(!is_clonal, Nmuts>0) %>% rename(n_sig=Nmuts)
      
      background_sigs = other_sigs %>%
        group_by(Signature) %>%
        # n total of mutations per cluster
        summarise(n_sig=sum(n_sig)) %>%
        mutate(n=sum(n_sig)) %>%
        mutate(Exposure=n_sig/n)

      tmp_final_table = lapply(unique(other_sigs$cluster_tool), function(c) {
        tmp = other_sigs %>% filter(cluster_tool == c)
        result_table = compare_signatures(df1=tmp, df2=clonal_sigs)
        result = result_table$df %>% mutate(match=result_table$match,
                                            cs_exp=result_table$cs_exp,
                                            cs_nsig=result_table$cs_n,
                                            n_diff=result_table$n_diff,
                                            n_tot=result_table$n_tot) %>%
          mutate(cluster_tool=as.character(c))
        return(result)
      }) %>% bind_rows()

      tmp_final_table_background = lapply(unique(other_sigs$cluster_tool), function(c){
        tmp = other_sigs %>% filter(cluster_tool == c)
        result_table = compare_signatures(df1=tmp, df2=background_sigs)
        result = result_table$df %>% mutate(match=result_table$match,
                                            cs_exp=result_table$cs_exp,
                                            cs_nsig=result_table$cs_n,
                                            n_diff=result_table$n_diff,
                                            n_tot=result_table$n_tot) %>%
          mutate(cluster_tool=as.character(c))
        return(result)
      }) %>% bind_rows()

      if (nrow(tmp_final_table) > 0) {
        f_background = tmp_final_table_background %>%
          mutate(n_rel=n_diff/n_tot) %>%
          select(cluster_tool, match, cs_exp, n_rel) %>%
          distinct() %>%
          dplyr::rename(bg_match=match, bg_cs_exp=cs_exp, bg_n_rel=n_rel) %>%
          filter(!is.na(cluster_tool))

        final_table_signature = tmp_final_table %>%
          mutate(n_rel=n_diff/n_tot) %>%
          select(cluster_tool, match, cs_exp, n_rel) %>%
          distinct() %>%
          left_join(f_background)
      }

      final_table_subclonal %>%
        left_join(final_table_signature) %>%
        # filter(nMuts > 100) %>%
        mutate(driver=ifelse(is_driver == F, 0, 1),
               bg_cs=1 - bg_cs_exp,
               cs_sign=1 - cs_exp,
               cs_sign=ifelse(is.na(cs_exp) & is_clonal == T, 1, cs_sign),
               bg_sign=ifelse(is.na(bg_cs) & is_clonal == T, 1, bg_cs),
               n_rel=ifelse(is.na(cs_exp) & is_clonal == T, 1, n_rel)) %>%
        rowwise() %>%
        mutate(
          score_all=(driver + n_never_tail + ((cs_sign+n_rel+bg_sign)/3))/3,
          score_driver=driver,
          score_sign=(cs_sign+n_rel+bg_sign)/3,
          score_tail=n_never_tail,
          score_no_tail=(driver + ((cs_sign+n_rel+bg_sign)/3))/2,
          score_no_driver=(((cs_sign+n_rel+bg_sign)/3) + n_never_tail)/2,
          score_no_sign=(driver + n_never_tail)/2) %>% 
        ungroup() %>% 
        mutate(tool=tool, signature_type=sign_type, patient_id=patient_id)
    }) %>% bind_rows()
  }) %>% bind_rows()
}) %>% bind_rows() %>% 
  select(patient_id, cluster_tool, everything())


pl_scores = score_table %>%
  pivot_longer(cols=c(score_driver, score_all, score_tail, score_no_driver, score_no_tail, score_no_sign, score_sign)) %>%
  mutate(name=factor(name, levels=c("score_driver", "score_tail", "score_sign","score_no_driver", "score_no_tail", "score_no_sign", "score_all"))) %>%
  ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0.9, ymax=1, fill="palegreen4", alpha=0.2) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0.55, ymax=.9, fill="goldenrod", alpha=0.2) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0.2, ymax=0.55, fill="salmon1", alpha=0.2) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0.2, ymax=0, fill="gainsboro", alpha=0.2) +
  geom_point(aes(x=name, y=value, col=cluster_tool), size=4) +
  geom_line(data=~filter(.x, !is.na(value)),
            aes(x=name, y=value, color=cluster_tool, group=cluster_tool), linewidth=.6)  +
  geom_text(data=~filter(.x, is_clonal & name=="score_all"),
            aes(x=name, y=value, color=cluster_tool, group=cluster_tool, label="Clonal"), vjust=0.07) +
  scale_color_manual("Cluster", values=colors_cluster) +
  scale_x_discrete(labels=c("score_driver"="Driver",
                            "score_tail"="Tail",
                            "score_sign"="Signature",
                            "score_no_driver"="Signature\nTail",
                            "score_no_tail"="Driver\nSignature",
                            "score_no_sign"="Driver\nTail",
                            "score_all"="All")) +
  facet_grid(patient_id ~ tool + signature_type) +
  theme_bw() + theme(axis.title.x=element_blank())

ggsave(pl_scores, filename=file.path(opt$prefix, "scores_clusters.pdf"), width=12, height=10)


pairs = combn(samples[samples != ''], 2, simplify = FALSE)
cluster_plots_pair = list()

for (i in 1:length(pairs)) {
  s1 = paste0('VAF.', pairs[[i]][1])
  s2 = paste0('VAF.', pairs[[i]][2])

  plot = data_viber %>%
    ggplot(aes(x =.data[[s1]], y=.data[[s2]], color=cluster)) +
    geom_point(alpha=0.2, size=.5) +
    scale_color_manual("Tool clusters", values = colors_cluster) +
    xlim(0,1) + ylim(0,1) +
    xlab(s1) + ylab(s2) +
    theme_bw() + ggtitle(tool) +
    guides(color=guide_legend(override.aes=list(size=3, alpha=1)))

  plot = plot + ggrepel::geom_label_repel(
    data=data_viber %>% filter(driver == TRUE),
    aes(x=.data[[s1]], y=.data[[s2]], label=gene, colour=cluster)),
    show.legend=F, inherit.aes=FALSE, size=3, min.segment.length=0,
    box.padding=1, max.overlaps=50)

  cluster_plots_pair[[i]] = plot
}

pl_multivariate = wrap_plots(cluster_plots_pair, ncol=3, guides="collect")

ggsave(pl_multivariate, filename=file.path(opt$prefix, "multivariate_clusters.pdf"), width=12, height=10)
